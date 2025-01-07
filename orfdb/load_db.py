#!/usr/bin/env python3
"""
Database loading script for the ORF database.

This script handles the loading of various annotation sources into the ORF database,
including GENCODE, RefSeq, CHESS, OpenProt, and custom ORF annotations.
"""

import logging
import time
import sys
import six
from os import listdir
from os.path import join, isfile
from collections import defaultdict
from pathlib import Path

import click
import pandas as pd

from orfdb import base, settings
from orfdb import annotation_loading, annotation_updates
from seqmap import genomic
import orfdb.util as orf_utils
import seqmap.utils as seq_utils
from sqlalchemy_batch_inserts import enable_batch_inserting


def configure_logger(log_file=None, level=logging.INFO, overwrite_log=True,
                    format=logging.BASIC_FORMAT):
    """Configure logging for the application.

    Args:
        log_file (str, optional): Path to log file. Defaults to None (stdout).
        level (int, optional): Logging level. Defaults to logging.INFO.
        overwrite_log (bool, optional): Whether to overwrite existing log.
            Defaults to True.
        format (str, optional): Log message format. Defaults to BASIC_FORMAT.
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
def load_db(drop_all):
    """Main entry point for loading the ORF database.

    Args:
        drop_all (bool): If True, drops and recreates all database tables.
    """
    configure_logger(
        f"{time.strftime('%Y%m%d_%H%M%S')}_orfdb_load_db.log",
        level=logging.INFO
    )

    if drop_all:
        logging.info('Dropping the database')
        base.Base.metadata.drop_all()
        
        logging.info('Creating the database')
        base.Base.metadata.create_all()

    # Initialize database session
    session = base.Session()
    enable_batch_inserting(session)

    logging.info(f'Loading {settings.db_connection_string}')

    logging.info('Loading GENCODE')
    #load_gencode_gff(session, settings.gencode_directory,
    #                 settings.gencode_version, settings.genomes_directory)
    
    logging.info('Loading Uniprot')
    #annotation_loading.load_uniprot(session, settings.uniprot_directory)
    
    logging.info('Loading RefSeq')
    #load_refseq_gff(session, settings.refseq_directory)

    logging.info('Loading CHESS')
    #load_chess_gff(session, settings.chess_directory)

    logging.info('Loading Openprot')
    #load_openprot_bed(session, settings.openprot_directory)
    
    logging.info('Updating RiboSeq CDS')
    #annotation_loading.update_riboseq_cds(session, settings.gencode_directory,
    #                                     'v42')

    logging.info('Loading Velia Phase 1 to 5 ORFs')
    #load_psl_phases_1to5(session, settings.velia_directory)

    logging.info('Loading Velia Phase 6 ORFs')
    #load_psl_phase_6(session, settings.velia_directory)

    logging.info('Updating Velia Phase ORFs')
    #update_psl_phases(session, settings.velia_directory) 

    logging.info('Updating ensembl <-> refseq gene mappings')
    #annotation_updates.update_ensembl_entrez_gene_mapping(
    #    settings.gencode_directory, 'v42', session)
    
    logging.info('Updating CHESS transcript IDs to fix parsing error')
    annotation_updates.update_chess_transcript_ids(
        settings.chess_directory, session)
    
    session.close()
    base.Session.close_all()


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
    gencode_expanded_gff = gencode_dir.joinpath(
        version,
        f'gencode.{version}.chr_patch_hapl_scaff.annotation.expanded.gff3'
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

    gff_df = pd.read_csv(gencode_expanded_gff, sep='\t', low_memory=False)
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
    session.upsert(
        base.Dataset,
        name='ENSEMBL',
        description='Automated ENSEMBL/GENCODE annotations',
        type='dataset',
        attrs={'version': version}
    )
    session.upsert(
        base.Dataset,
        name='HAVANA',
        description='Manual ENSEMBL/GENCODE annotations',
        type='dataset',
        attrs={'version': version}
    )
    session.upsert(
        base.Dataset,
        name='CCDS',
        description='CCDS consensus cds project',
        type='dataset'
    )
    session.upsert(
        base.Dataset,
        name='HGNC_ID',
        description='HUGO Human gene nomenclature numeric gene IDs',
        type='dataset',
        attrs={'version': version}
    )
    session.upsert(
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


def load_openprot_bed(session, openprot_dir):
    """Load OpenProt annotations from BED files.

    Args:
        session: SQLAlchemy session object
        openprot_dir (Path): Directory containing OpenProt files
    """
    openprot_bed_phased = openprot_dir.joinpath(
        'v1.6',
        'human-openprot-r1_6-refprots+altprots+isoforms-_phased.bed'
    )
    
    bed_df = pd.read_csv(openprot_bed_phased, sep='\t', index_col=0)

    # Create dataset
    openprot_dataset = session.upsert(
        base.Dataset,
        name="openprot",
        description="Openprot database",
        type="dataset",
        attrs={"version": "v1.6"}
    )

    # Load OpenProt components
    logging.info('Loading Openprot CDS')
    annotation_loading.load_openprot_cds(session, bed_df)

    logging.info('Loading Openprot ORFs')
    cds_orf_map = annotation_loading.load_openprot_orfs(session, bed_df)
    
    logging.info('Loading Openprot CDS <-> ORF mappings')
    annotation_loading.load_refseq_orf_cds(session, cds_orf_map)


def load_psl_phases_1to5(session, velia_dir):
    """Load Phase 1-5 PSL data.

    Args:
        session: SQLAlchemy session object
        velia_dir (Path): Directory containing Phase 1-5 files
    """
    psl_df, phase_df = orf_utils.etl_phase1to5_psl(session, velia_dir)

    logging.info('Loading Velia CDS from PSL')
    annotation_loading.load_psl_phase_cds(session, psl_df, 'velia')

    logging.info('Loading Velia phase ORFs from PSL')
    annotation_loading.load_psl_phase_orfs_legacy(
        session, psl_df, phase_df, 'velia')
    
    logging.info('Loading Velia phase proteins from phase df')
    annotation_loading.load_psl_phase_proteins(session, phase_df, 'velia')


def load_psl_phase_6(session, velia_dir):
    """Load Phase 6 PSL data.

    Args:
        session: SQLAlchemy session object
        velia_dir (Path): Directory containing Phase 6 files
    """
    psl_df, phase_df = orf_utils.etl_phase6_psl(session, velia_dir)

    logging.info('Loading Velia CDS from PSL')
    annotation_loading.load_psl_phase_cds(session, psl_df, 'velia')

    logging.info('Loading Velia phase ORFs from PSL')
    annotation_loading.load_psl_phase_orfs_legacy(
        session, psl_df, phase_df, 'velia')

    logging.info('Loading Velia phase proteins from phase df')
    annotation_loading.load_psl_phase_proteins(session, phase_df, 'velia')


def update_psl_phases(session, velia_dir):
    """Update PSL phase data.

    Args:
        session: SQLAlchemy session object
        velia_dir (Path): Directory containing phase data files
    """
    psl_phase1to5_df, _ = orf_utils.etl_phase1to5_psl(session, velia_dir)
    psl_phase6_df, _ = orf_utils.etl_phase6_psl(session, velia_dir)

    logging.info('Updating Velia phase 1 to 5 ORFs from PSL')
    annotation_loading.update_psl_phase_orfs(session, psl_phase1to5_df)

    logging.info('Updating Velia phase 6 ORFs from PSL')
    annotation_loading.update_psl_phase_orfs(session, psl_phase6_df)


if __name__ == '__main__':
    load_db()
