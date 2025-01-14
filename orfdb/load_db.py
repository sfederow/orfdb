#!/usr/bin/env python3
"""
Database loading script for the ORF database.

This script handles the loading of various annotation sources into the ORF database,
including GENCODE, RefSeq, CHESS, OpenProt, and custom ORF annotations.
"""

import logging
import time
import sys
from pathlib import Path
from typing import Optional

import click
import pandas as pd
from sqlalchemy import select
from sqlalchemy.orm import Session

from orfdb import base, settings
from orfdb import annotation_loading, annotation_updates
from seqmap import genomic
import orfdb.util as orf_utils
import seqmap.utils as seq_utils
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

    # Initialize database session
    session = base.Session()
    enable_batch_inserting(session)

    logging.info(f'Loading {settings.db_connection_string}')

    try:
        logging.info('Loading GENCODE')
        load_gencode_gff(session, settings.gencode_directory,
                        settings.gencode_version, settings.genomes_directory)

        logging.info('Loading RefSeq')
        load_refseq_gff(session, settings.refseq_directory)

        logging.info('Loading CHESS')
        load_chess_gff(session, settings.chess_directory)

        logging.info('Loading BigProt')
        load_bigprot_tables(session, settings.bigprot_directory)

        logging.info('Updating ensembl <-> refseq gene mappings')
        annotation_updates.update_chess_transcript_ids(
            settings.chess_directory, session)
        
        session.commit()
    except Exception as e:
        session.rollback()
        logging.error(f"Error during database loading: {str(e)}")
        raise
    finally:
        session.close()


def load_gencode_gff(session, gencode_dir, version, genome_dir):
    """Load GENCODE annotations from GFF files.

    Args:
        session: SQLAlchemy session object
        gencode_dir (Path): Directory containing GENCODE files
        version (str): GENCODE version (e.g., 'v42')
        genome_dir (Path): Directory containing genome files
    """
    assembly_info = gencode_dir.joinpath(
        version,
        'GCA_000001405.28_GRCh38.p13_assembly_report.txt'
    )
    gencode_gff = gencode_dir.joinpath(
        version,
        f'gencode.{version}.chr_patch_hapl_scaff.annotation.gff3.gz'
    )

    assembly_df = pd.read_csv(
        assembly_info,
        sep='\t',
        comment='#',
        names=['Sequence-Name', 'Sequence-Role', 'Assigned-Molecule',
               'Assigned-Molecule-Location/Type', 'GenBank-Accn',
               'Relationship', 'RefSeq-Accn', 'Assembly-Unit',
               'Sequence-Length', 'UCSC-style-name']
    )

    logging.info(f'Loading GENCODE {version} GFF')

    # Read gzipped GFF file, skipping header rows
    gff_df = pd.read_csv(
        gencode_gff, 
        sep='\t', 
        low_memory=False,
        compression='gzip' if str(gencode_gff).endswith('.gz') else None,
        comment='#',
        names=[
            'seqid', 'source', 'type', 'start', 'end', 
            'score', 'strand', 'phase', 'attributes'
        ]
    )
    gff_df.fillna('', inplace=True)

    # Filter dataframes by type
    gene_gff_df = gff_df[gff_df['type'] == 'gene'].copy()
    tx_gff_df = gff_df[gff_df['type'] == 'transcript'].copy()
    exon_gff_df = gff_df[gff_df['type'] == 'exon'].copy()
    cds_gff_df = gff_df[gff_df['type'] == 'CDS'].copy()
    utr_five_gff_df = gff_df[gff_df['type'].isin(['five_prime_UTR'])].copy()
    utr_three_gff_df = gff_df[gff_df['type'].isin(['three_prime_UTR'])].copy()

    del gff_df

    # Create or update datasets
    base.upsert(
        session,
        base.Dataset,
        name='ENSEMBL',
        description='Automated ENSEMBL/GENCODE annotations',
        type='dataset',
        attrs={'version': version}
    )
    base.upsert(
        session,
        base.Dataset,
        name='HAVANA',
        description='Manual ENSEMBL/GENCODE annotations',
        type='dataset',
        attrs={'version': version}
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
        attrs={'version': version}
    )
    base.upsert(
        session,
        base.Dataset,
        name='HGNC_SYMBOL',
        description='HUGO Human gene nomenclature symbols',
        type='dataset',
        attrs={'version': version}
    )

    logging.info('Loading genome assembly')
    annotation_loading.load_genome_assembly(
        session, assembly_df, assembly_info.stem)

    # Map assembly IDs
    assembly_ids = {}
    for assembly in session.query(base.Assembly).all():
        if assembly.ucsc_style_name.startswith('chr') and \
           len(assembly.ucsc_style_name) < 6:
            assembly_ids[assembly.ucsc_style_name] = assembly.id
        else:
            assembly_ids[assembly.genbank_accession] = assembly.id

    # Load various GENCODE components
    logging.info('Loading GENCODE genes')
    annotation_loading.load_gencode_genes(session, gene_gff_df, assembly_ids)

    logging.info('Loading GENCODE exons')
    annotation_loading.load_gencode_exons(session, exon_gff_df, assembly_ids)

    logging.info('Loading GENCODE transcripts')
    transcript_exon_map = annotation_loading.load_gencode_transcripts(
        session, tx_gff_df, exon_gff_df, gencode_dir, version, assembly_ids)

    logging.info('Loading GENCODE transcript <-> exon relationships')
    annotation_loading.load_gencode_transcript_exons(
        session, transcript_exon_map)

    logging.info('Loading GENCODE UTRs')
    annotation_loading.load_gencode_utr(
        session, utr_five_gff_df, assembly_ids, base.UtrFive)
    annotation_loading.load_gencode_utr(
        session, utr_three_gff_df, assembly_ids, base.UtrThree)

    logging.info('Loading GENCODE cds')
    annotation_loading.load_gencode_cds(session, cds_gff_df, assembly_ids)
    
    logging.info('Loading GENCODE ORFs')
    cds_orf_map = annotation_loading.load_gencode_orfs(
        session, cds_gff_df, assembly_ids)
    annotation_loading.load_gencode_orf_cds(session, cds_orf_map)

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
        gencode_dir.joinpath(version, 'Ribo-seq_ORFs.bed'),
        names=bed_cols,
        sep='\t'
    )

    genome_seq = genomic.load_genome(
        genome_dir.joinpath('hg38', 'hg38_no-altChr.fa'))

    logging.info('Loading GENCODE Riboseq ORFs')
    (transcript_exon_rb_map,
     transcript_orf_map,
     exon_cds_map,
     cds_orf_map) = annotation_loading.load_gencode_riboseq_orfs(
        session, orf_bed_df, assembly_ids, genome_seq)

    annotation_loading.load_gencode_riboseq_orf_cds(session, cds_orf_map)

    logging.info('Loading GENCODE transcript <-> exon relationships for '
                'riboseq orfs')
    annotation_loading.load_gencode_riboseq_transcript_exons(
        session, transcript_exon_rb_map)
    
    logging.info('Loading GENCODE transcript <-> ORF relationships')
    annotation_loading.load_gencode_riboseq_transcript_orfs(
        session, transcript_orf_map)
    annotation_loading.load_gencode_transcript_orfs(
        session, tx_gff_df, assembly_ids)

    logging.info('Loading GENCODE exon <-> cds relationships')
    annotation_loading.load_gencode_riboseq_exon_cds(
        session, exon_cds_map)
    annotation_loading.load_gencode_exon_cds(
        session, cds_gff_df, assembly_ids)

    # Clean up
    del gene_gff_df
    del tx_gff_df
    del exon_gff_df
    del cds_gff_df
    del utr_five_gff_df
    del utr_three_gff_df


def load_refseq_gff(session, refseq_dir):
    """Load RefSeq annotations from GFF files.

    Args:
        session: SQLAlchemy session object
        refseq_dir (Path): Directory containing RefSeq files
    """
    refseq_gff = refseq_dir.joinpath(
        'GCA_000001405.15_GRCh38_full_analysis_set.'
        'refseq_annotation.expanded.ids.gff'
    )
    refseq_df = pd.read_csv(
        refseq_gff,
        sep='\t',
        dtype={'HGNC_ID': str, 'entrez_gene_id': str}
    )
    refseq_df.fillna('', inplace=True)

    refseq_cols = [
        'seq_id', 'source', 'type', 'start', 'end', 'score', 'strand',
        'phase', 'ID', 'Dbxref', 'Name', 'chromosome', 'assembly_id',
        'HGNC_ID', 'entrez_gene_id', 'gbkey', 'genome', 'mol_type',
        'description', 'gene', 'gene_biotype', 'Parent', 'product',
        'transcript_id', 'gene_synonym', 'model_evidence', 'tag',
        'protein_id', 'experiment', 'inference', 'Note'
    ]

    refseq_df = refseq_df[refseq_cols].copy()

    # Filter dataframes by type
    gene_gff_df = refseq_df[refseq_df['type'] == 'gene'].copy()
    tx_gff_df = refseq_df[refseq_df['type'].isin(
        ['mRNA', 'lnc_RNA', 'transcript'])].copy()
    exon_gff_df = refseq_df[refseq_df['type'] == 'exon'].copy()
    cds_gff_df = refseq_df[refseq_df['type'] == 'CDS'].copy()

    # Create dataset
    refseq_dataset = session.upsert(
        base.Dataset,
        name="RefSeq",
        description="NCBI Genome Annotation from gff",
        type="dataset",
        attrs={"version": "GCA_000001405.15_GRCh38"}
    )

    for source in set(refseq_df['source']):
        session.upsert(
            base.Dataset,
            name=source,
            description='',
            type="dataset",
            attrs={"version": "GCA_000001405.15_GRCh38"}
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

    logging.info('Loading RefSeq ORFs')
    cds_orf_map = annotation_loading.load_refseq_orfs(session, cds_gff_df)

    logging.info('Loading RefSeq ORFs <-> CDS mapping')
    annotation_loading.load_refseq_orf_cds(session, cds_orf_map)

    logging.info('Loading RefSeq transcript <-> orf mappings')
    annotation_loading.load_refseq_transcript_orfs(session, cds_gff_df)


def load_chess_gff(session, chess_dir):
    """Load CHESS annotations from GFF files.

    Args:
        session: SQLAlchemy session object
        chess_dir (Path): Directory containing CHESS files
    """
    chess_expanded_gff = chess_dir.joinpath('chess3.0.expanded.gff')
    chess_df = pd.read_csv(chess_expanded_gff, sep='\t', low_memory=False)
    chess_df.fillna('', inplace=True)

    # Create dataset
    chess_dataset = session.upsert(
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
        session.upsert(
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


def load_bigprot_tables(session, bigprot_dir):
    """Load BigProt tables from the given directory.

    Args:
        session: SQLAlchemy session object
        bigprot_dir (Path): Directory containing BigProt files
    """

    required_files = [
        'orfset_v0.9.3_full_minlen_15_maxlen_999999999_orfs.csv.gz',
        'orfset_v0.9.3_full_minlen_15_maxlen_999999999_transcript_orfs.csv.gz',
        'orfset_v0.9.3_full_minlen_15_maxlen_999999999_cds_orfs.csv.gz'
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
                verbose=True,
                dataset_name='BigProt',
                genome_fasta_fpath=str(settings.genomes_directory.joinpath(
                    'hg38', 'GRCh38.p13.genome.fa.gz')),
                db_settings_fpath=str(Path(__file__).parent.parent / 'settings.ini'),
                gtf_file_fpath=str(bigprot_dir.joinpath(
                    'gencode.v42.chr_patch_hapl_scaff.annotation.expanded.gff3')),
                min_codon_length=15,
                max_codon_length=999999999
            )
        except Exception as e:
            logging.error(f'Failed to run BigProt analysis: {str(e)}')
            return
    
    try:
        logging.info('Loading BigProt ORFs')
        annotation_loading.load_bigprot_orfs(session, bigprot_dir)
        
        logging.info('Loading BigProt transcripts')
        annotation_loading.load_bigprot_transcripts(session, bigprot_dir)
        
        logging.info('Loading BigProt CDS-ORF mappings')
        annotation_loading.load_bigprot_cds_orf(session, bigprot_dir)
        
        session.commit()
    except Exception as e:
        session.rollback()
        logging.error(f'Failed to load BigProt data: {str(e)}')
        raise


if __name__ == '__main__':
    try:
        load_db()
    except Exception as e:
        logging.error(f"Failed to load database: {str(e)}")
        sys.exit(1)
