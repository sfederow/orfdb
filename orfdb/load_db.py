#!/usr/bin/env python3
"""
Database loading script for the ORF database.

This script handles the loading of various annotation sources into the ORF database,
including GENCODE, RefSeq, CHESS, OpenProt, and custom ORF annotations.
"""

import logging, re, sys, time
from pathlib import Path
from typing import Optional

import click
import pandas as pd
from sqlalchemy import select
from sqlalchemy.orm import Session

from orfdb import base, settings
from orfdb import annotation_loading, annotation_updates
import orfdb.util as orf_utils
from sqlalchemy_batch_inserts import enable_batch_inserting


def configure_logger(
    log_file: Optional[str] = None, 
    level: int = logging.INFO, 
    overwrite_log: bool = True,
    format: str = logging.BASIC_FORMAT
) -> None:
    """Configure logging for the application.

    Args:
        log_file: Path to log file. Defaults to None (stdout).
        level: Logging level. Defaults to logging.INFO.
        overwrite_log: Whether to overwrite existing log. Defaults to True.
        format: Log message format. Defaults to BASIC_FORMAT.
    """
    if log_file is None:
        logging.basicConfig(stream=sys.stdout, level=level, format=format)
    else:
        logging.basicConfig(
            filename=log_file,
            level=level,
            filemode=('w' if overwrite_log else 'a'),
            format=format
        )
        console = logging.StreamHandler()
        console.setLevel(level)
        console.setFormatter(logging.Formatter(format))
        logging.getLogger('').addHandler(console)


def bump_version(version: str) -> str:
    """
    Bumps the patch version by 1 for a semantic version string (e.g., v0.9.3 -> v0.9.4).
    
    Args:
        version (str): The version string to bump (e.g., 'v0.9.3').
    
    Returns:
        str: The bumped version string.
    """
    
    if version is None:
        return 'v0.0.1'

    match = re.match(r"^(v?)(\d+)\.(\d+)\.(\d+)$", version)
    if not match:
        raise ValueError(f"Invalid version format: {version}")
    
    prefix, major, minor, patch = match.groups()
    # Increment the patch number
    bumped_version = f"{prefix}{major}.{minor}.{int(patch) + 1}"
    return bumped_version


@click.command()
@click.option('--drop-all', is_flag=True, help='Empty database and reload data.')
def load_db(drop_all: bool) -> None:
    """Main entry point for loading the ORF database.

    Args:
        drop_all: If True, drops and recreates all database tables.
    """
    configure_logger(
        f"{time.strftime('%Y%m%d_%H%M%S')}_orfdb_load_db.log",
        level=logging.INFO
    )

    if drop_all:
        logging.info('Dropping the database')
        base.Base.metadata.drop_all(bind=base.engine)
        
        logging.info('Creating the database')
        base.Base.metadata.create_all(bind=base.engine)

    session = base.Session()
    enable_batch_inserting(session)

    logging.info(f'Loading {settings.db_connection_string}')

    try:
        logging.info('Loading genome assembly')
        load_genome_assembly(session, settings.data_dir, settings.genome_assembly)

        logging.info('Loading GENCODE')
        load_gencode_gff(session, settings.data_dir, settings.gencode_gff, settings.gencode_refseq)

        logging.info('Loading RefSeq')
        load_refseq_gff(session, settings.data_dir, settings.refseq_gff)

        logging.info('Loading CHESS')
        load_chess_gff(session, settings.data_dir, settings.chess_gff)

        logging.info('Loading BigProt')
        #load_bigprot_tables(session, settings.data_dir, settings.genome, settings.bigprot_version, new_run=False)
        
        session.commit()

    except Exception as e:
        session.rollback()
        logging.error(f"Error during database loading: {str(e)}")
        raise
    finally:
        session.close()


def load_genome_assembly(session, data_dir, genome_assembly_file):
    """
    Load genome assembly information into the database.

    This function reads the genome assembly report from a specified directory,
    processes the data, and loads  it into the database using the provided session.

    Args:
        session: SQLAlchemy session object used for database transactions.
        genome_dir (Path): Directory containing the genome assembly report file.

    The function expects the genome assembly report to be located at:
    'hg38/GCF_000001405.40.assembly_report.tsv' within the provided genome_dir.
    """
    assembly_info = data_dir.joinpath(genome_assembly_file)

    column_names = [
        'Sequence-Name', 
        'Sequence-Role', 
        'Assigned-Molecule',
        'Assigned-Molecule-Location/Type', 
        'GenBank-Accn', 
        'Relationship', 
        'RefSeq-Accn', 
        'Assembly-Unit', 
        'Sequence-Length', 
        'UCSC-style-name'
    ]

    assembly_df = pd.read_csv(
        assembly_info,
        sep='\t',
        comment='#',
        names=column_names
    )

    annotation_loading.load_genome_assembly(
        session, assembly_df, assembly_info.stem)


def load_gencode_gff(session, data_dir, gencode_gff_file, gencode_refseq_file):
    """Load GENCODE annotations from GFF files.

    Args:
        session: SQLAlchemy session object
        gencode_dir (Path): Directory containing GENCODE files
    """


    logging.info(f'Loading GENCODE {gencode_gff_file} GFF')

    gff_df = pd.read_csv(
        data_dir.joinpath(gencode_gff_file), 
        sep='\t', 
        low_memory=False,
        compression='gzip' if str(gencode_gff_file).endswith('.gz') else None,
        comment='#',
        names=[
            'seqid', 'source', 'type', 'start', 'end', 
            'score', 'strand', 'phase', 'attributes'
        ]
    )
    gff_df.fillna('', inplace=True)


    assembly_ids = {}
    for assembly in session.query(base.Assembly).all():
        if assembly.ucsc_style_name.startswith('chr') and \
           len(assembly.ucsc_style_name) < 6:
            assembly_ids[assembly.ucsc_style_name] = assembly.id
        else:
            assembly_ids[assembly.genbank_accession] = assembly.id

    gff_df['assembly_id'] = gff_df.apply(
        lambda x: assembly_ids[x.seqid] if x.seqid in assembly_ids.keys() else '', 
        axis=1
    )

    gff_df = gff_df[gff_df['assembly_id'] != ''].copy()

    gene_gff_df = gff_df[gff_df['type'] == 'gene'].copy()
    tx_gff_df = gff_df[gff_df['type'] == 'transcript'].copy()
    exon_gff_df = gff_df[gff_df['type'] == 'exon'].copy()
    cds_gff_df = gff_df[gff_df['type'] == 'CDS'].copy()
    utr_five_gff_df = gff_df[gff_df['type'].isin(['five_prime_UTR'])].copy()
    utr_three_gff_df = gff_df[gff_df['type'].isin(['three_prime_UTR'])].copy()

    del gff_df

    base.upsert(
        session,
        base.Dataset,
        name='ENSEMBL',
        description='Automated ENSEMBL/GENCODE annotations',
        type='dataset',
        attrs={'version': gencode_gff_file.stem}
    )
    base.upsert(
        session,
        base.Dataset,
        name='HAVANA',
        description='Manual ENSEMBL/GENCODE annotations',
        type='dataset',
        attrs={'version': gencode_gff_file.stem}
    )
    base.upsert(
        session,
        base.Dataset,
        name='CCDS',
        description='CCDS consensus cds project',
        type='dataset'
    )
    base.upsert(
        session,
        base.Dataset,
        name='HGNC_ID',
        description='HUGO Human gene nomenclature numeric gene IDs',
        type='dataset',
        attrs={'version': gencode_gff_file.stem}
    )
    base.upsert(
        session,
        base.Dataset,
        name='HGNC_SYMBOL',
        description='HUGO Human gene nomenclature symbols',
        type='dataset',
        attrs={'version': gencode_gff_file.stem}
    )

    logging.info('Loading GENCODE genes')
    annotation_loading.load_gencode_genes(session, gene_gff_df)

    logging.info('Loading GENCODE exons')
    annotation_loading.load_gencode_exons(session, exon_gff_df)

    logging.info('Loading GENCODE transcripts')
    transcript_exon_map = annotation_loading.load_gencode_transcripts(
        session, tx_gff_df, exon_gff_df, data_dir.joinpath(gencode_refseq_file))

    logging.info('Loading GENCODE transcript <-> exon relationships')
    annotation_loading.load_gencode_transcript_exons(
        session, transcript_exon_map)

    logging.info('Loading GENCODE UTRs')
    annotation_loading.load_gencode_utr(
        session, utr_five_gff_df, base.UtrFive)
    annotation_loading.load_gencode_utr(
        session, utr_three_gff_df, base.UtrThree)

    logging.info('Loading GENCODE cds')
    annotation_loading.load_gencode_cds(session, cds_gff_df)
    
    logging.info('Loading GENCODE exon <-> cds relationships')
    annotation_loading.load_gencode_exon_cds(
        session, cds_gff_df)

    # Free up memory
    del gene_gff_df
    del tx_gff_df
    del exon_gff_df
    del cds_gff_df
    del utr_five_gff_df
    del utr_three_gff_df


def load_refseq_gff(session, data_dir, refseq_gff_file):
    """Load RefSeq annotations from GFF files.

    Args:
        session: SQLAlchemy session object
        refseq_dir (Path): Directory containing RefSeq files
    """

    refseq_version = '.'.join(str(refseq_gff_file).split('.')[0:2])

    refseq_df = pd.read_csv(
        data_dir.joinpath(refseq_gff_file),
        sep='\t',
        dtype={'HGNC_ID': str, 'entrez_gene_id': str},
        compression='gzip' if str(refseq_gff_file).endswith('.gz') else None,
        comment='#',
        names=[
            'seq_id', 'source', 'type', 'start', 'end', 'score', 'strand',
            'phase', 'attributes'
        ]
    )
    refseq_df.fillna('', inplace=True)

    assembly_ids = {}
    for assembly in session.query(base.Assembly).all():
        assembly_ids[assembly.refseq_accession] = assembly.id

    refseq_df['assembly_id'] = refseq_df['seq_id'].apply(
        lambda x: int(assembly_ids.get(x, -1))
    )

    refseq_df = refseq_df[refseq_df['assembly_id'] != -1].copy()

    # Filter dataframes by type
    gene_gff_df = refseq_df[refseq_df['type'] == 'gene'].copy()
    tx_gff_df = refseq_df[refseq_df['type'].isin(
        ['mRNA', 'lnc_RNA', 'transcript'])].copy()
    exon_gff_df = refseq_df[refseq_df['type'] == 'exon'].copy()
    cds_gff_df = refseq_df[refseq_df['type'] == 'CDS'].copy()

    # Create dataset
    refseq_dataset = base.upsert(
        session,
        base.Dataset,
        name="RefSeq",
        description="NCBI Genome Annotation from gff",
        type="dataset",
        attrs={"version": refseq_version}
    )

    for source in set(refseq_df['source']):
        base.upsert(
            session,
            base.Dataset,
            name=source,
            description='',
            type="dataset",
            attrs={"version": refseq_version}
        )

    # Load RefSeq components
    logging.info('Loading RefSeq genes')
    annotation_loading.load_refseq_genes(session, gene_gff_df)
    
    logging.info('Loading RefSeq exons')
    annotation_loading.load_refseq_exons(session, exon_gff_df)

    logging.info('Loading RefSeq transcripts')
    transcript_exon_map = annotation_loading.load_refseq_transcripts(
        session, tx_gff_df, exon_gff_df)

    logging.info('Loading RefSeq transcript <-> exon mappings')
    annotation_loading.load_refseq_transcript_exons(
        session, transcript_exon_map)

    logging.info('Loading RefSeq CDS')
    annotation_loading.load_refseq_cds(session, cds_gff_df)


def load_chess_gff(session, data_dir, chess_gff_file):
    """Load CHESS annotations from GFF files.

    Args:
        session: SQLAlchemy session object
        chess_dir (Path): Directory containing CHESS files
    """

    chess_df = pd.read_csv(
        data_dir.joinpath(chess_gff_file),
        sep='\t',
        compression='gzip' if str(chess_gff_file).endswith('.gz') else None,
        low_memory=False,
        comment='#',
        names=[
            'seq_id', 'source', 'type', 'start', 'end', 'score', 'strand',
            'phase', 'attributes'
        ]
    )
    chess_df.fillna('', inplace=True)

    assembly_ids = {}
    for assembly in session.query(base.Assembly).all():
        if assembly.ucsc_style_name.startswith('chr') and \
           len(assembly.ucsc_style_name) < 6:
            assembly_ids[assembly.ucsc_style_name] = assembly.id
        else:
            assembly_ids[assembly.genbank_accession] = assembly.id

    chess_df['assembly_id'] = chess_df['seq_id'].apply(
        lambda x: int(assembly_ids.get(x, -1))
    )

    chess_df = chess_df[chess_df['assembly_id'] != -1].copy()

    # Create dataset
    chess_dataset = base.upsert(
        session,
        base.Dataset,
        name="CHESS",
        description="CHESS annotation source",
        type="dataset",
        attrs={
            "version": "GCA_000001405.15_GRCh38",
            "chess_version": 3.0
        }
    )

    for source in set(chess_df['source']):
        base.upsert(
            session,
            base.Dataset,
            name=source,
            description='',
            type="dataset",
            attrs={
                "version": "GCA_000001405.15_GRCh38",
                "chess_version": 3.0
            }
        )
    
    # Filter dataframes by type
    tx_gff_df = chess_df[chess_df['type'] == 'transcript'].copy()
    exon_gff_df = chess_df[chess_df['type'] == 'exon'].copy()
    cds_gff_df = chess_df[chess_df['type'] == 'CDS'].copy()

    # Load CHESS components
    logging.info('Loading CHESS exons')
    annotation_loading.load_chess_exons(session, exon_gff_df)
    
    logging.info('Loading CHESS transcripts')
    transcript_exon_map = annotation_loading.load_chess_transcripts(
        session, tx_gff_df, exon_gff_df)

    logging.info('Loading CHESS transcrips <-> exons')
    annotation_loading.load_chess_transcript_exons(
        session, transcript_exon_map)

    logging.info('Loading CHESS CDS')
    annotation_loading.load_chess_cds(session, cds_gff_df)


def load_bigprot_tables(session, data_dir, genome_file, bigprot_version=None, new_run=False):
    """Load BigProt tables from the given directory.

    Args:
        session: SQLAlchemy session object
        bigprot_dir (Path): Directory containing BigProt files
    """

    if not new_run and bigprot_version:
        return

    bigprot_version = bump_version(bigprot_version)
    bigprot_dir = data_dir.joinpath('bigprot', bigprot_version)

    required_files = [
        f'orfset_BigProt_minlen_15_maxlen_999999999_orfs.csv.gz',
        f'orfset_BigProt_minlen_15_maxlen_999999999_transcript_orfs.csv.gz',
        f'orfset_BigProt_minlen_15_maxlen_999999999_cds_orfs.csv.gz'
    ]

    files_exist = all(
        bigprot_dir.joinpath(filename).exists() 
        for filename in required_files
    )

    if not files_exist:
        logging.info('BigProt analysis files missing, running analysis...')
        try:
            from orfdb.bigprot.find_orfs import perform_analysis
            perform_analysis(
                output=str(bigprot_dir),
                verbose=False,
                dataset_name='BigProt',
                genome_fasta_fpath=str(data_dir.joinpath(genome_file)),
                db_settings_fpath=str(Path(__file__).parent.parent / 'settings.ini'),
                min_codon_length=15,
                max_codon_length=999999999,
                skip_annotation=True,
                transcript_chunk_size=1000,
                max_chunks=1000,
                num_processes=10,
                accession_namespace='genbank',
            )
        except Exception as e:
            logging.error(f'Failed to run BigProt analysis: {str(e)}')
            return
    
    try:
        logging.info('Loading BigProt ORFs')
        annotation_loading.load_bigprot_orfs(session, bigprot_dir)
        
        logging.info('Loading BigProt transcripts')
        annotation_loading.load_bigprot_transcript_orfs(session, bigprot_dir)
        
        logging.info('Loading BigProt CDS-ORF mappings')
        annotation_loading.load_bigprot_cds_orf(session, bigprot_dir)
        
        session.commit()
    except Exception as e:
        session.rollback()
        logging.error(f'Failed to load BigProt data: {str(e)}')
        raise


def load_gencode_riboseq(session, gencode_riboseq_file):
    """Load GENCODE Ribo-seq ORFs from a file.

    Args:
        session: SQLAlchemy session object
        gencode_riboseq_file (Path): Path to the GENCODE Ribo-seq ORFs file
    """

    gencode_riboseq_file = settings.data_dir.joinpath(gencode_riboseq_file)

    logging.info(f'Loading GENCODE Ribo-seq ORFs from {gencode_riboseq_file}')

    # Load Ribo-seq ORFs
    bed_cols = [
        'chrom', 'chromStart', 'chromEnd', 'name', 'score',
        'strand', 'thickStart', 'thickEnd', 'itemRgb',
        'blockCount', 'blockSizes', 'chromStarts', 'name2',
        'cdsStartStat', 'cdsEndStat', 'exonFrames', 'type',
        'geneName', 'geneName2', 'geneType', 'transcript_biotype',
        'sequence', 'all_transcript_ids', 'all_gene_ids',
        'replicated', 'ref_studies'
    ]

    orf_bed_df = pd.read_csv(
        settings.data_dir.joinpath(gencode_riboseq_file),
        names=bed_cols,
        sep='\t'
    )

    logging.info('Load annotation for GENCODE Riboseq ORFs')
    annotation_loading.load_gencode_riboseq(session, orf_bed_df)



if __name__ == '__main__':
    try:
        load_db()
    except Exception as e:
        logging.error(f"Failed to load database: {str(e)}")
        sys.exit(1)
