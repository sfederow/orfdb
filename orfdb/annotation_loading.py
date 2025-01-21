"""Module to load core annotation data.

This module provides functions for loading various types of annotation data
into the ORF database, including genome assemblies, genes, transcripts,
exons, UTRs, and CDS regions from different data sources.
"""

from Bio.Seq import Seq
from orfdb.base import *  
from orfdb import base
from orfdb import util as orf_utils  

import ast
import gzip
import hashlib
import logging
import pandas as pd
from os.path import basename
from warnings import warn
from pathlib import Path
from sqlalchemy.orm import Session
from sqlalchemy import select, or_, and_
from typing import Dict, List, Tuple



def load_genome_assembly(session: Session, assembly_df: pd.DataFrame, assembly_name: str) -> None:
    """Load genome assembly information into the database.

    Args:
        session: SQLAlchemy session
        assembly_df: DataFrame containing assembly information
        assembly_name: Name of the assembly
    """
    # Add unmapped assembly if it doesn't exist
    stmt = select(base.Assembly).filter_by(genbank_accession='unmapped')
    unmapped = session.execute(stmt).scalar_one_or_none()
    
    if not unmapped:
        unmapped = base.Assembly(
            genbank_accession='unmapped',
            refseq_accession='unmapped',
            ucsc_style_name='unmapped',
            sequence_role='unmapped',
            assembly_unit='unmapped',
            assigned_molecule='unmapped',
            assigned_molecule_location='unmapped',
            sequence_length=-1,
            genome_accession='unmapped'
        )
        session.add(unmapped)
        session.commit()

    # Process each assembly row
    assemblies_to_add = []
    for _, row in assembly_df.iterrows():
        # Check if assembly exists
        stmt = select(base.Assembly).filter(
            or_(
                base.Assembly.genbank_accession == row['GenBank-Accn'],
                base.Assembly.refseq_accession == row['RefSeq-Accn']
            )
        )
        existing = session.execute(stmt).scalar_one_or_none()

        if existing:
            # Update existing assembly
            existing.genbank_accession = row['GenBank-Accn']
            existing.refseq_accession = row['RefSeq-Accn']
            existing.ucsc_style_name = row['UCSC-style-name']
            existing.sequence_role = row['Sequence-Role']
            existing.assembly_unit = row['Assembly-Unit']
            existing.assigned_molecule = row['Assigned-Molecule']
            existing.assigned_molecule_location = row['Assigned-Molecule-Location/Type']
            existing.sequence_length = int(row['Sequence-Length'])
            existing.genome_accession = assembly_name
            existing.attrs = {
                'sequence_name': row['Sequence-Name'],
                'relationship': row['Relationship']
            }
        else:
            # Create new assembly
            assembly = base.Assembly(
                genbank_accession=row['GenBank-Accn'],
                refseq_accession=row['RefSeq-Accn'],
                ucsc_style_name=row['UCSC-style-name'],
                sequence_role=row['Sequence-Role'],
                assembly_unit=row['Assembly-Unit'],
                assigned_molecule=row['Assigned-Molecule'],
                assigned_molecule_location=row['Assigned-Molecule-Location/Type'],
                sequence_length=int(row['Sequence-Length']),
                genome_accession=assembly_name,
                attrs={
                    'sequence_name': row['Sequence-Name'],
                    'relationship': row['Relationship']
                }
            )
            assemblies_to_add.append(assembly)

    # Bulk insert new assemblies
    if assemblies_to_add:
        session.bulk_save_objects(assemblies_to_add)
        session.commit()

    logging.info(f'Added {len(assemblies_to_add)} genome assemblies')


def load_gencode_genes(
    session: Session, 
    gene_gff_df: pd.DataFrame, 
    assembly_ids: Dict[str, int]
) -> None:
    """Load GENCODE gene annotations into the database.

    Args:
        session: SQLAlchemy session
        gene_gff_df: DataFrame containing gene annotations
        assembly_ids: Dictionary mapping chromosome names to assembly IDs
    """
    # Map assembly IDs to DataFrame
    gene_gff_df['assembly_id'] = gene_gff_df.apply(
        lambda x: assembly_ids[x.seqid], 
        axis=1
    )

    # Group by unique columns to handle redundant entries
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    gene_gff_df.sort_values(by=unique_cols, inplace=True)
    
    # Process attributes before grouping
    gene_gff_df['attrs'] = gene_gff_df['attributes'].apply(parse_attributes)
    gene_gff_df['hgnc_id'] = gene_gff_df['attrs'].apply(lambda x: x.get('hgnc_id', ''))
    gene_gff_df['gene_name'] = gene_gff_df['attrs'].apply(lambda x: x.get('gene_name', ''))
    gene_gff_df['gene_id'] = gene_gff_df['attrs'].apply(lambda x: x.get('gene_id', ''))
    gene_gff_df['gene_type'] = gene_gff_df['attrs'].apply(lambda x: x.get('gene_type', ''))
    
    grouped_gene_gff_df = gene_gff_df.groupby(unique_cols).agg({
        'hgnc_id': list,
        'gene_name': list,
        'gene_id': list,
        'gene_type': 'first',
        'source': list,
        'attrs': list
    })

    # Prepare genes for bulk insert
    genes_to_add = []
    synonym_dict = {}

    for idx, row in grouped_gene_gff_df.iterrows():
        gene = base.Gene(
            start=idx[0],
            end=idx[1],
            strand=idx[2],
            assembly_id=idx[3],
            hgnc_id=row.hgnc_id[0],
            hgnc_name=row.gene_name[0],
            ensembl_id=row.gene_id[0],
            refseq_id='',
            chess_id='',
            velia_id='',
            gene_type=row.gene_type,
            long_name='',
            attrs={}
        )
        genes_to_add.append(gene)
        
        # Build synonym dictionary for xrefs
        syn_dict = {source: set() for source in row.source}
        for i, gene_id in enumerate(row.gene_id):
            syn_dict[row.source[i]].add(gene_id)
        
        syn_dict['HGNC_ID'] = set(row.hgnc_id)
        syn_dict['HGNC_SYMBOL'] = set(row.gene_name)
        synonym_dict[idx] = syn_dict

    # Add unmapped gene
    stmt = select(base.Assembly).filter_by(genbank_accession='unmapped')
    unmapped_assembly = session.execute(stmt).scalar_one()
    
    genes_to_add.append(base.Gene(
        start=-1,
        end=-1,
        strand='',
        assembly_id=unmapped_assembly.id,
        hgnc_id="unmapped",
        hgnc_name="unmapped",
        ensembl_id="unmapped",
        refseq_id="unmapped",
        chess_id="unmapped",
        velia_id="unmapped",
        gene_type="unmapped",
        long_name="",
        attrs={}
    ))

    # Bulk insert genes
    session.bulk_save_objects(genes_to_add)
    session.commit()
    
    logging.info(f'Added {len(genes_to_add)} GENCODE genes')

    # Get gene IDs for xrefs
    stmt = select(base.Gene)

    gene_map = {
        (g.start, g.end, g.strand, g.assembly_id): g.id 
        for g in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs
    stmt = select(base.Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }

    # Prepare xrefs for bulk insert
    xrefs_to_add = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                xrefs_to_add.append(base.SequenceRegionXref(
                    sequence_region_id=gene_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['ENSEMBL'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    # Bulk insert xrefs
    session.bulk_save_objects(xrefs_to_add)
    session.commit()
    
    logging.info(f'Added {len(xrefs_to_add)} GENCODE gene synonyms')


def load_gencode_exons(session: Session, exon_gff_df: pd.DataFrame, assembly_ids: Dict[str, int]) -> None:
    """Load exons from GENCODE GFF file.

    Args:
        session: SQLAlchemy session object
        exon_gff_df: GENCODE GFF dataframe filtered for exons
        assembly_ids: Mapping of sequence IDs to assembly IDs
    """
    # Map assembly IDs to DataFrame
    exon_gff_df['assembly_id'] = exon_gff_df.apply(
        lambda x: assembly_ids[x.seqid],
        axis=1
    )

    # Group by unique columns to handle redundant entries
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    exon_gff_df.sort_values(by=unique_cols, inplace=True)
    
    # Process attributes before grouping
    exon_gff_df['attrs'] = exon_gff_df['attributes'].apply(parse_attributes)
    exon_gff_df['exon_id'] = exon_gff_df['attrs'].apply(lambda x: x.get('exon_id', ''))
    exon_gff_df['ID'] = exon_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    
    grouped_exon_gff_df = exon_gff_df.groupby(unique_cols).agg({
        'exon_id': list,
        'ID': list,
        'source': list,
        'attrs': list
    })
        
    exons = []
    synonym_dict = {}

    for idx, row in grouped_exon_gff_df.iterrows():
        exons.append(base.Exon(
            start=idx[0],
            end=idx[1],
            strand=idx[2],
            assembly_id=idx[3],
            ensembl_id=row.exon_id[0],
            refseq_id='',
            chess_id='',
            velia_id=''
        ))

        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.exon_id[i])
            syn_dict[row.source[i]].add(ID)

        synonym_dict[idx] = syn_dict

    # Bulk insert exons
    session.bulk_save_objects(exons)
    session.commit()

    logging.info(f'Added {len(exons)} GENCODE exons')
    
    # Get exon IDs for xrefs
    stmt = select(base.Exon)

    exon_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id 
        for e in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs
    stmt = select(base.Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }
    
    # Prepare xrefs for bulk insert
    exon_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                exon_xrefs.append(base.SequenceRegionXref(
                    sequence_region_id=exon_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['ENSEMBL'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    # Bulk insert xrefs
    session.bulk_save_objects(exon_xrefs)
    session.commit()

    logging.info(f'Added {len(exon_xrefs)} GENCODE exon synonyms')


def load_gencode_transcripts(
    session: Session, 
    tx_gff_df: pd.DataFrame, 
    exon_gff_df: pd.DataFrame, 
    gencode_dir: Path,
    version: str, 
    assembly_ids: Dict[str, int]
) -> Dict[str, List[Tuple[int, int]]]:
    """Load transcripts from GENCODE GFF file.

    Args:
        session: SQLAlchemy session object
        tx_gff_df: GENCODE GFF dataframe filtered for transcripts
        exon_gff_df: GENCODE GFF dataframe filtered for exons
        gencode_dir: Path to GENCODE directory
        version: GENCODE version
        assembly_ids: Mapping of sequence IDs to assembly IDs

    Returns:
        Dictionary mapping transcript IDs to exon information
    """
    # Map assembly IDs
    tx_gff_df['assembly_id'] = tx_gff_df.apply(
        lambda x: assembly_ids[x.seqid],
        axis=1
    )
    exon_gff_df['assembly_id'] = exon_gff_df.apply(
        lambda x: assembly_ids[x.seqid],
        axis=1
    )

    tx_gff_df['attrs'] = tx_gff_df['attributes'].apply(parse_attributes)
    tx_gff_df['ID'] = tx_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    tx_gff_df['gene_id'] = tx_gff_df['attrs'].apply(lambda x: x.get('gene_id', ''))
    tx_gff_df['havana_transcript'] = tx_gff_df['attrs'].apply(lambda x: x.get('havana_transcript', ''))
    tx_gff_df['ccdsid'] = tx_gff_df['attrs'].apply(lambda x: x.get('ccdsid', ''))
    tx_gff_df['transcript_support_level'] = tx_gff_df['attrs'].apply(lambda x: x.get('transcript_support_level', ''))
    tx_gff_df['transcript_type'] = tx_gff_df['attrs'].apply(lambda x: x.get('transcript_type', ''))
    tx_gff_df['tag'] = tx_gff_df['attrs'].apply(lambda x: x.get('tag', ''))

    exon_gff_df['attrs'] = exon_gff_df['attributes'].apply(parse_attributes)
    exon_gff_df['ID'] = exon_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    exon_gff_df['transcript_id'] = exon_gff_df['attrs'].apply(lambda x: x.get('transcript_id', ''))
    exon_gff_df['exon_number'] = exon_gff_df['attrs'].apply(lambda x: x.get('exon_number', ''))


    # Merge transcript and exon data
    merged_df = tx_gff_df.merge(
        exon_gff_df,
        left_on=('ID', 'seqid'),
        right_on=('transcript_id', 'seqid'),
        suffixes=('_tx', '_ex'),
        how='left'
    )

    # Group by transcript information
    transcript_exon_gff_df = merged_df.groupby(
        ['ID_tx', 'start_tx', 'end_tx', 'strand_tx', 'assembly_id_tx']
    ).aggregate(list)

    # Get existing exon and gene mappings using SQLAlchemy 2.0 style
    stmt = select(base.Exon)
    ensembl_orf_exon_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id
        for e in session.execute(stmt).scalars().all()
    }

    stmt = select(base.Gene)
    ensembl_orf_gene_map = {
        g.ensembl_id: g.id
        for g in session.execute(stmt).scalars().all()
    }

    # Load GENCODE-RefSeq mapping
    gencode_refseq_map = {}
    with gzip.open(
        gencode_dir.joinpath(
            version,
            f'gencode.{version}.metadata.RefSeq.gz'
        ), 'rt'
    ) as infile:
        for line in infile.readlines():
            vals = line.rstrip('\n').split('\t')
            gencode_refseq_map[vals[0]] = vals[1]

    transcripts = []
    transcript_exon_map = {}
    synonym_dict = {}
    existing_entries = []

    # Get unmapped gene using SQLAlchemy 2.0 style
    stmt = select(base.Gene).filter_by(hgnc_id='unmapped')
    unmapped_gene = session.execute(stmt).scalar_one()

    for (transcript_id, start, end, strand, assembly_id), row in \
            transcript_exon_gff_df.iterrows():
        block_sizes = []
        chrom_starts = []
        exon_ids = []

        try:
            gene_id = ensembl_orf_gene_map[row.gene_id[0]]
        except KeyError:
            gene_id = unmapped_gene.id

        try:
            refseq_id = gencode_refseq_map[transcript_id]
        except KeyError:
            refseq_id = ''

        for i in range(len(row.ID_ex)):
            if pd.isna(row.ID_ex[i]):
                continue
            chrom_starts.append(str(int(row.start_ex[i])))
            block_sizes.append(str(int(row.end_ex[i] - row.start_ex[i])))
            exon_ids.append((
                row.exon_number[i],
                ensembl_orf_exon_map[
                    (row.start_ex[i],
                     row.end_ex[i],
                     row.strand_ex[i],
                     assembly_id)
                ]
            ))

        block_sizes = ';'.join(block_sizes)
        chrom_starts = ';'.join(chrom_starts)
        transcript_idx_str = (
            f'{start}_{end}_{strand}_{int(assembly_id)}_'
            f'{block_sizes}_{chrom_starts}'
        )
        transcript_idx = hashlib.sha3_512(
            transcript_idx_str.encode('utf-8')
        ).hexdigest()
        
        if transcript_idx in existing_entries:
            synonym_dict[transcript_idx]['ENSEMBL'].update([transcript_id])
            synonym_dict[transcript_idx]['HAVANA'].update(
                row.havana_transcript
            )
            synonym_dict[transcript_idx]['CCDS'].update(row.ccdsid)
        else:
            existing_entries.append(transcript_idx)
            synonym_dict[transcript_idx] = {}
            synonym_dict[transcript_idx]['ENSEMBL'] = set([transcript_id])
            synonym_dict[transcript_idx]['HAVANA'] = set(
                row.havana_transcript
            )
            synonym_dict[transcript_idx]['CCDS'] = set(row.ccdsid)

            transcripts.append(base.Transcript(
                start=start,
                end=end,
                strand=strand,
                assembly_id=assembly_id,
                block_sizes=block_sizes,
                chrom_starts=chrom_starts,
                transcript_idx=transcript_idx,
                transcript_idx_str=transcript_idx_str,
                gene_id=gene_id,
                ensembl_id=transcript_id,
                refseq_id=refseq_id,
                chess_id='',
                velia_id='',
                support_level=row.transcript_support_level[0],
                transcript_type=row.transcript_type[0],
                attrs={'tag': row.tag[0]}
            ))

            transcript_exon_map[transcript_idx] = []

        for exon_num, exon_id in exon_ids:
            transcript_exon_map[transcript_idx].append((exon_num, exon_id))

    # Bulk insert transcripts
    session.bulk_save_objects(transcripts)
    session.commit()

    logging.info(f'Added {len(transcripts)} GENCODE transcripts')

    # Get transcript IDs for xrefs using SQLAlchemy 2.0 style
    stmt = select(base.Transcript)
    transcript_map = {
        t.transcript_idx: t.id
        for t in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs using SQLAlchemy 2.0 style
    stmt = select(base.Dataset)
    dataset_ids = {
        d.name: d.id
        for d in session.execute(stmt).scalars().all()
    }

    # Prepare xrefs for bulk insert
    transcript_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                transcript_xrefs.append(base.SequenceRegionXref(
                    sequence_region_id=transcript_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['ENSEMBL'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    # Bulk insert xrefs
    session.bulk_save_objects(transcript_xrefs)
    session.commit()

    logging.info(f'Added {len(transcript_xrefs)} GENCODE transcript synonyms')

    return transcript_exon_map


def load_gencode_utr(
    session: Session, 
    utr_gff_df: pd.DataFrame, 
    assembly_ids: Dict[str, int], 
    utr_class
) -> None:
    """Load UTRs from GENCODE GFF file.

    Args:
        session: SQLAlchemy session object
        utr_gff_df: Subset of the GENCODE GFF with UTR records
        assembly_ids: Mapping of GENCODE 'seqid' to internal database IDs
        utr_class: The UTR class to use (either UTR3 or UTR5)
    """
    # Map assembly IDs to DataFrame
    utr_gff_df['assembly_id'] = utr_gff_df.apply(
        lambda x: assembly_ids[x.seqid], 
        axis=1
    )

    # Group by unique columns to handle redundant entries
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    utr_gff_df.sort_values(by=unique_cols, inplace=True)
    
    # Process attributes before grouping
    utr_gff_df['attrs'] = utr_gff_df['attributes'].apply(parse_attributes)
    utr_gff_df['ID'] = utr_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    
    grouped_utr_gff_df = utr_gff_df.groupby(unique_cols).agg({
        'ID': list,
        'source': list,
        'attrs': list
    })
        
    utrs = []
    synonym_dict = {}

    for idx, row in grouped_utr_gff_df.iterrows():
        utrs.append(utr_class(
            start=idx[0],
            end=idx[1],
            strand=idx[2],
            assembly_id=idx[3],
            ensembl_id=row.ID[0],
            refseq_id='',
            chess_id='',
            velia_id=''
        ))
        
        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(ID)

        synonym_dict[idx] = syn_dict

    # Bulk insert UTRs
    session.bulk_save_objects(utrs)
    session.commit()

    logging.info(f'Added {len(utrs)} GENCODE UTRs')
    
    # Get UTR IDs for xrefs
    stmt = select(utr_class)
    utr_map = {
        (u.start, u.end, u.strand, u.assembly_id): u.id 
        for u in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs
    stmt = select(base.Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }
    
    # Prepare xrefs for bulk insert
    utr_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                utr_xrefs.append(base.SequenceRegionXref(
                    sequence_region_id=utr_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['ENSEMBL'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    # Bulk insert xrefs
    session.bulk_save_objects(utr_xrefs)
    session.commit()

    logging.info(f'Added {len(utr_xrefs)} GENCODE UTR synonyms')


def load_gencode_cds(
    session: Session, 
    cds_gff_df: pd.DataFrame, 
    assembly_ids: Dict[str, int]
) -> None:
    """Load CDS from GENCODE GFF file.

    Args:
        session: SQLAlchemy session object
        cds_gff_df: GENCODE GFF dataframe filtered for CDS records
        assembly_ids: Mapping of sequence IDs to assembly IDs
    """
    # Map assembly IDs to DataFrame
    cds_gff_df['assembly_id'] = cds_gff_df.apply(
        lambda x: assembly_ids[x.seqid], 
        axis=1
    )

    # Group by unique columns to handle redundant entries
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    cds_gff_df.sort_values(by=unique_cols, inplace=True)
    
    # Process attributes before grouping
    cds_gff_df['attrs'] = cds_gff_df['attributes'].apply(parse_attributes)
    cds_gff_df['ID'] = cds_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    cds_gff_df['protein_id'] = cds_gff_df['attrs'].apply(lambda x: x.get('protein_id', ''))
    cds_gff_df['Parent'] = cds_gff_df['attrs'].apply(lambda x: x.get('Parent', ''))
    
    grouped_cds_gff_df = cds_gff_df.groupby(unique_cols).agg({
        'ID': list,
        'protein_id': list,
        'Parent': list,
        'source': list,
        'attrs': list
    })

    stmt = select(base.Gene).filter_by(hgnc_id='unmapped')
    unmapped_gene = session.execute(stmt).scalar_one()

    # Get transcript mappings
    stmt = select(base.Transcript)
    transcript_map = {
        t.ensembl_id: t.id 
        for t in session.execute(stmt).scalars().all()
    }

    cds_entries = []
    synonym_dict = {}

    for idx, row in grouped_cds_gff_df.iterrows():
        try:
            transcript_id = transcript_map[row.Parent[0]]
        except KeyError:
            transcript_id = unmapped_gene.id

        cds_entries.append(base.Cds(
            start=idx[0],
            end=idx[1],
            strand=idx[2],
            assembly_id=idx[3],
            ensembl_id=row.ID[0],
            ensembl_protein_id=row.protein_id[0] if row.protein_id[0] else '',
            refseq_id='',
            chess_id='',
            velia_id='',
        ))


        syn_dict = {source: set() for source in row.source}
        
        for i, ID in enumerate(row.ID):
            if ID:
                syn_dict[row.source[i]].add(ID)
            if row.protein_id[i]:
                syn_dict[row.source[i]].add(row.protein_id[i])

        synonym_dict[idx] = syn_dict

    # Bulk insert CDS entries
    session.bulk_save_objects(cds_entries)
    session.commit()

    logging.info(f'Added {len(cds_entries)} GENCODE CDS')

    # Get CDS IDs for xrefs
    stmt = select(base.Cds)
    cds_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id 
        for c in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs
    stmt = select(base.Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }

    cds_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                cds_xrefs.append(base.SequenceRegionXref(
                    sequence_region_id=cds_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['ENSEMBL'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    session.bulk_save_objects(cds_xrefs)
    session.commit()

    logging.info(f'Added {len(cds_xrefs)} GENCODE CDS synonyms')


def load_gencode_exon_cds(session: Session, cds_gff_df: pd.DataFrame, assembly_ids: Dict[str, int]) -> None:
    """Load exon-CDS relationships from GENCODE data.

    Args:
        session: SQLAlchemy session object
        cds_gff_df: GENCODE GFF dataframe filtered for CDS records
        assembly_ids: Mapping of sequence IDs to assembly IDs
    """

    stmt = select(base.Cds)
    ensembl_vdb_cds_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id 
        for c in session.execute(stmt).scalars().all()
    }

    stmt = select(base.SequenceRegionXref, base.Exon).join(
        base.Exon, 
        base.Exon.id == base.SequenceRegionXref.sequence_region_id
    )
    ensembl_vdb_exon_map = {
        (exon.assembly_id, syn.xref): exon.id 
        for syn, exon in session.execute(stmt).all()
    }

    cds_gff_df['assembly_id'] = cds_gff_df.apply(lambda x: assembly_ids[x.seqid], axis=1)
    cds_gff_df['attrs'] = cds_gff_df['attributes'].apply(parse_attributes)
    cds_gff_df['exon_id'] = cds_gff_df['attrs'].apply(lambda x: x.get('exon_id', ''))

    exon_cds_ids = []
    exon_cds = []

    for i, row in cds_gff_df.iterrows():
        exon_cds_ids.append((
            ensembl_vdb_exon_map[(row.assembly_id, row.exon_id)], 
            ensembl_vdb_cds_map[(row.start, row.end, row.strand, row.assembly_id)]
        ))

    stmt = select(base.ExonCds)
    existing = {(ec.exon_id, ec.cds_id) for ec in session.execute(stmt).scalars().all()}
    
    new_relationships = set(exon_cds_ids).difference(existing)

    for (exon_id, cds_id) in new_relationships:
        exon_cds.append(base.ExonCds(
            exon_id=exon_id,
            cds_id=cds_id
        ))

    session.add_all(exon_cds)
    session.commit()

    logging.info(f'Added {len(exon_cds)} GENCODE exon-CDS relationships')


def load_gencode_transcript_exons(
    session: Session, 
    transcript_exon_map: Dict[str, List[Tuple[int, int]]]
) -> None:
    """Load transcript-exon relationships from GENCODE data.

    Args:
        session: SQLAlchemy session object
        transcript_exon_map: Mapping of transcript IDs to exon info (number, id)
    """

    stmt = select(base.Transcript)
    transcript_map = {
        t.transcript_idx: t.id
        for t in session.execute(stmt).scalars().all()
    }

    transcript_exons = []
    seen_relationships = set()

    for transcript_idx, exon_entries in transcript_exon_map.items():
        transcript_id = transcript_map[transcript_idx]
        
        for exon_num, exon_id in exon_entries:
            # Create unique key to avoid duplicates
            relationship_key = (transcript_id, exon_id, exon_num)
            if relationship_key in seen_relationships:
                continue
                
            seen_relationships.add(relationship_key)
            
            transcript_exons.append(base.TranscriptExon(
                transcript_id=transcript_id,
                exon_id=exon_id,
                exon_number=exon_num
            ))

    if transcript_exons:
        session.bulk_save_objects(transcript_exons)
        session.commit()

    logging.info(f'Added {len(transcript_exons)} GENCODE transcript-exon relationships')


def load_uniprot(session, uniprot_dir):
    """Load UniProt protein data.

    Args:
        session: SQLAlchemy session object
        uniprot_dir (Path): Directory containing UniProt files
    """
    uniprot_df = pd.read_csv(
        uniprot_dir.joinpath('uniprot_sprot_human.tab'),
        sep='\t'
    )

    # Create dataset
    uniprot_dataset = session.upsert(
        Dataset,
        name="uniprot",
        description="UniProt protein database",
        type="dataset",
        attrs={"version": "2023_02"}
    )

    proteins = []
    for i, row in uniprot_df.iterrows():
        proteins.append(Protein(
            row.Entry,
            row.Entry_name,
            row.Status,
            row.Protein_names,
            row.Gene_names,
            row.Organism,
            row.Length,
            row.Sequence,
            row.Gene_names_primary,
            row.EC_number,
            row.Gene_names_synonym,
            row.Gene_names_ordered_locus,
            row.Gene_names_ORF,
            row.Protein_names_intermediate,
            row.Protein_names_short,
            row.Protein_names_EC,
            row.Protein_names_allergen,
            row.Protein_names_CD_antigen,
            row.Protein_names_INN,
            row.Proteomes,
            row.Proteomes_component
        ))

    session.add_all(proteins)
    session.commit()

    logging.info(f'Added {len(proteins)} UniProt proteins')


def load_refseq_genes(session: Session, gene_gff_df: pd.DataFrame) -> None:
    """Load genes from RefSeq GFF file.

    Args:
        session: SQLAlchemy session object
        gene_gff_df: RefSeq GFF dataframe filtered for gene records
    """

    gene_gff_df['attrs'] = gene_gff_df['attributes'].apply(parse_attributes)
    gene_gff_df['ID'] = gene_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    gene_gff_df['Name'] = gene_gff_df['attrs'].apply(lambda x: x.get('Name', ''))
    gene_gff_df['description'] = gene_gff_df['attrs'].apply(lambda x: x.get('description', ''))
    gene_gff_df['gene_biotype'] = gene_gff_df['attrs'].apply(lambda x: x.get('gene_biotype', ''))
    gene_gff_df['HGNC_ID'] = gene_gff_df['attrs'].apply(lambda x: x.get('HGNC_ID', ''))
    gene_gff_df['entrez_gene_id'] = gene_gff_df['attrs'].apply(lambda x: x.get('entrez_gene_id', ''))

    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    gene_gff_df.sort_values(by=unique_cols, inplace=True)

    grouped_gene_gff_df = gene_gff_df.groupby(unique_cols).aggregate(list)

    stmt = select(base.Gene)
    vdb_gene_idx_map = {
        (g.start, g.end, g.strand, g.assembly_id): g.id 
        for g in session.execute(stmt).scalars().all()
    }

    update_entries = []
    refseq_only_genes = []
    synonym_dict = {}

    for gene_idx, row in grouped_gene_gff_df.iterrows():
        if gene_idx in vdb_gene_idx_map.keys():
            gene_id = vdb_gene_idx_map[gene_idx]

            update_entries.append({
                "id": gene_id, 
                "refseq_id": row.entrez_gene_id[0],
                "attrs": {
                    "description": row.description[0],
                    "gene_biotype": row.gene_biotype[0]
                }
            })

        elif gene_idx not in synonym_dict.keys():

            refseq_only_genes.append(base.Gene(
                start=gene_idx[0], 
                end=gene_idx[1], 
                strand=gene_idx[2], 
                assembly_id=gene_idx[3],
                hgnc_id=row.HGNC_ID[0],
                hgnc_name=row.Name[0],
                ensembl_id='',
                refseq_id=row.entrez_gene_id[0],
                chess_id='',
                velia_id='',
                gene_type=row.gene_biotype[0],
                attrs={'description': row.description[0]}
            ))
                    
        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.Name[i])
            syn_dict[row.source[i]].add(ID)
            syn_dict[row.source[i]].add(row.entrez_gene_id[i])

        syn_dict['HGNC_ID'] = set(row.HGNC_ID)
        synonym_dict[gene_idx] = syn_dict

    if update_entries:
        session.bulk_update_mappings(base.Gene, update_entries)
        session.commit()

    if refseq_only_genes:
        session.bulk_save_objects(refseq_only_genes)
        session.commit()

    logging.info(f'Updated {len(update_entries)} genes with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(refseq_only_genes)} RefSeq genes without an exact coordinate match to GENCODE')

    stmt = select(base.Gene)
    vdb_gene_idx_map = {
        (g.start, g.end, g.strand, g.assembly_id): g.id 
        for g in session.execute(stmt).scalars().all()
    }    

    stmt = select(base.Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }

    gene_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                gene_xrefs.append(base.SequenceRegionXref(
                    sequence_region_id=vdb_gene_idx_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['RefSeq'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    if gene_xrefs:
        session.bulk_save_objects(gene_xrefs)
        session.commit()

    logging.info(f'Added {len(gene_xrefs)} RefSeq gene synonyms')


def load_refseq_transcripts(session: Session, tx_gff_df: pd.DataFrame, exon_gff_df: pd.DataFrame) -> Dict[str, List[Tuple[int, int]]]:
    """Load transcripts from RefSeq GFF file.

    Args:
        session: SQLAlchemy session object
        tx_gff_df: RefSeq GFF dataframe filtered for transcript records
        exon_gff_df: RefSeq GFF dataframe filtered for exon records

    Returns:
        Dict[str, List[Tuple[int, int]]]: Mapping of transcript IDs to lists of (exon_number, exon_id) tuples
    """

    tx_gff_df['attrs'] = tx_gff_df['attributes'].apply(parse_attributes)
    tx_gff_df['ID'] = tx_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    tx_gff_df['gene'] = tx_gff_df['attrs'].apply(lambda x: x.get('gene', ''))
    tx_gff_df['HGNC_ID'] = tx_gff_df['attrs'].apply(lambda x: x.get('HGNC_ID', ''))
    tx_gff_df['entrez_gene_id'] = tx_gff_df['attrs'].apply(lambda x: x.get('entrez_gene_id', ''))
    tx_gff_df['product'] = tx_gff_df['attrs'].apply(lambda x: x.get('product', ''))
    tx_gff_df['model_evidence'] = tx_gff_df['attrs'].apply(lambda x: x.get('model_evidence', ''))
    tx_gff_df['type'] = tx_gff_df['attrs'].apply(lambda x: x.get('type', ''))
    tx_gff_df['Parent'] = tx_gff_df['attrs'].apply(lambda x: x.get('Parent', ''))
    tx_gff_df['inference'] = tx_gff_df['attrs'].apply(lambda x: x.get('inference', ''))
    tx_gff_df['Note'] = tx_gff_df['attrs'].apply(lambda x: x.get('Note', ''))
    tx_gff_df['experiment'] = tx_gff_df['attrs'].apply(lambda x: x.get('experiment', ''))
    tx_gff_df['gbkey'] = tx_gff_df['attrs'].apply(lambda x: x.get('gbkey', ''))

    exon_gff_df['attrs'] = exon_gff_df['attributes'].apply(parse_attributes)
    exon_gff_df['ID'] = exon_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    exon_gff_df['gene'] = exon_gff_df['attrs'].apply(lambda x: x.get('gene', ''))
    exon_gff_df['HGNC_ID'] = exon_gff_df['attrs'].apply(lambda x: x.get('HGNC_ID', ''))
    exon_gff_df['entrez_gene_id'] = exon_gff_df['attrs'].apply(lambda x: x.get('entrez_gene_id', ''))
    exon_gff_df['product'] = exon_gff_df['attrs'].apply(lambda x: x.get('product', ''))
    exon_gff_df['model_evidence'] = exon_gff_df['attrs'].apply(lambda x: x.get('model_evidence', ''))
    exon_gff_df['type'] = exon_gff_df['attrs'].apply(lambda x: x.get('type', ''))
    exon_gff_df['Parent'] = exon_gff_df['attrs'].apply(lambda x: x.get('Parent', ''))
    exon_gff_df['inference'] = exon_gff_df['attrs'].apply(lambda x: x.get('inference', ''))
    exon_gff_df['Note'] = exon_gff_df['attrs'].apply(lambda x: x.get('Note', ''))
    exon_gff_df['experiment'] = exon_gff_df['attrs'].apply(lambda x: x.get('experiment', ''))
    exon_gff_df['gbkey'] = exon_gff_df['attrs'].apply(lambda x: x.get('gbkey', ''))
    
    exon_numbers = []
    for i, row in exon_gff_df.iterrows():
        vals = row.ID.split('-')
        if len(vals) < 3:
            exon_numbers.append(1)
        else:
            try:
                exon_numbers.append(int(vals[-1]))
            except:
                exon_numbers.append(1)
                
    exon_gff_df['exon_number'] = exon_numbers

    map_cols = ['seq_id', 'assembly_id', 'start', 'end', 'strand', 
                'gene', 'HGNC_ID', 'ID', 'entrez_gene_id', 
                'product', 'model_evidence', 'type', 'Parent',  
                'inference', 'Note', 'experiment', 'gbkey'] 

    merged_df = tx_gff_df[map_cols].merge(
        exon_gff_df[map_cols + ['exon_number']], 
        left_on=('ID', 'seq_id'), 
        right_on=('Parent', 'seq_id'),
        suffixes=('_tx', '_ex'), 
        how='left'
    )    

    merged_df.sort_values(by=['ID_tx', 'start_ex'], inplace=True)

    transcript_exon_gff_df = merged_df.groupby(['ID_tx', 'start_tx', 'end_tx', 'strand_tx', 'assembly_id_tx']).aggregate(list)

    stmt = select(Exon)
    vdb_exon_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id 
        for e in session.execute(stmt).scalars().all()
    }

    stmt = select(Gene)
    refseq_vdb_gene_map = {
        g.refseq_id: g.id 
        for g in session.execute(stmt).scalars().all()
    }

    stmt = select(Transcript)
    vdb_transcript_idx_map = {
        t.transcript_idx: t.id 
        for t in session.execute(stmt).scalars().all()
    }

    transcripts = []
    transcript_exon_map = {}
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    stmt = select(Gene).filter_by(hgnc_id='unmapped')
    unmapped_gene = session.execute(stmt).scalar_one()

    for (transcript_id, start, end, strand, assembly_id), row in transcript_exon_gff_df.iterrows():
        block_sizes = []
        chrom_starts = []
        exon_ids = []
        assembly_id = int(assembly_id)

        try:
            gene_id = refseq_vdb_gene_map[row.entrez_gene_id_tx[0]]
        except:
            gene_id = unmapped_gene.id

        for i in range(len(row.ID_ex)):
            if pd.isna(row.ID_ex[i]):
                continue
                
            chrom_starts.append(str(row.start_ex[i]))
            block_sizes.append(str(row.end_ex[i] - row.start_ex[i]))
            exon_ids.append((row.exon_number[i], vdb_exon_map[(row.start_ex[i], row.end_ex[i], row.strand_ex[i], assembly_id)]))

        block_sizes = ';'.join(block_sizes)
        chrom_starts = ';'.join(chrom_starts)
        transcript_idx_str = f'{start}_{end}_{strand}_{assembly_id}_{block_sizes}_{chrom_starts}'
        transcript_idx = hashlib.sha3_512(transcript_idx_str.encode('utf-8')).hexdigest()
        
        if transcript_idx in vdb_transcript_idx_map.keys():
            transcript_vdb_id = vdb_transcript_idx_map[transcript_idx]
            
            update_entries.append({
                "id": transcript_vdb_id, 
                "refseq_id": transcript_id,
                "attrs": {
                    "inference": row.inference_tx[0],
                    "experiment": row.experiment_tx[0],
                    "product": row['product_tx'][0],
                    "parent": row.Parent_tx[0],
                    "gbkey": row.gbkey_tx[0]
                }
            })

        elif transcript_idx not in existing_entries:
            existing_entries.append(transcript_idx)
            
            transcripts.append(Transcript(
                start=start, 
                end=end, 
                strand=strand, 
                assembly_id=assembly_id,
                block_sizes=block_sizes, 
                chrom_starts=chrom_starts, 
                transcript_idx=transcript_idx, 
                transcript_idx_str=transcript_idx_str,
                gene_id=gene_id, 
                ensembl_id='', 
                refseq_id=transcript_id, 
                chess_id='', 
                velia_id='',
                support_level=row.model_evidence_tx[0], 
                transcript_type=row.type_tx[0],
                attrs={
                    "inference": row.inference_tx[0],
                    "experiment": row.experiment_tx[0],
                    "product": row['product_tx'][0],
                    "parent": row.Parent_tx[0],
                    "gbkey": row.gbkey_tx[0]
                }
            ))
        
        if transcript_idx in synonym_dict.keys():
            synonym_dict[transcript_idx]['RefSeq'].update([transcript_id, transcript_id[4:]])
            synonym_dict[transcript_idx]['RefSeq'].update(row.Parent_ex)
            synonym_dict[transcript_idx]['HGNC_ID'].update(row.HGNC_ID_tx)
            synonym_dict[transcript_idx]['HGNC_SYMBOL'].update(row.gene_tx)
        else:
            synonym_dict[transcript_idx] = {}
            synonym_dict[transcript_idx]['RefSeq'] = set([transcript_id, transcript_id[4:]])
            synonym_dict[transcript_idx]['RefSeq'].update(row.Parent_ex)
            synonym_dict[transcript_idx]['HGNC_ID'] = set(row.HGNC_ID_tx)
            synonym_dict[transcript_idx]['HGNC_SYMBOL'] = set(row.gene_tx)

        transcript_exon_map[transcript_id] = []
        for (exon_num, exon_id) in exon_ids:
            transcript_exon_map[transcript_id].append((exon_num, exon_id))

    if update_entries:
        session.bulk_update_mappings(Transcript, update_entries)
        session.commit()

    if transcripts:
        session.bulk_save_objects(transcripts)
        session.commit()

    logging.info(f'Updated {len(update_entries)} transcripts with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(transcripts)} RefSeq transcripts without an exact coordinate match to GENCODE')

    stmt = select(Transcript)
    ensembl_vdb_tx_map = {
        t.transcript_idx: t.id 
        for t in session.execute(stmt).scalars().all()
    }

    stmt = select(Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }

    transcript_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                transcript_xrefs.append(TranscriptXref(
                    transcript_id=ensembl_vdb_tx_map[idx],
                    xref=synonym,
                    type='synonym',
                    transcript_dataset_id=dataset_ids['RefSeq'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    if transcript_xrefs:
        session.bulk_save_objects(transcript_xrefs)
        session.commit()

    logging.info(f'Added {len(transcript_xrefs)} RefSeq transcript synonyms')

    return transcript_exon_map


def load_refseq_exons(session: Session, exon_gff_df: pd.DataFrame) -> None:
    """Load exons from RefSeq GFF file.

    Args:
        session: SQLAlchemy session object
        exon_gff_df: RefSeq GFF dataframe filtered for exon records
    """
    # Group by unique columns to handle redundant entries
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    exon_gff_df.sort_values(by=unique_cols, inplace=True)

    # Process attributes
    exon_gff_df['attrs'] = exon_gff_df['attributes'].apply(parse_attributes)
    exon_gff_df['ID'] = exon_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    exon_gff_df['Name'] = exon_gff_df['attrs'].apply(lambda x: x.get('Name', ''))
    exon_gff_df['Parent'] = exon_gff_df['attrs'].apply(lambda x: x.get('Parent', ''))
    exon_gff_df['Note'] = exon_gff_df['attrs'].apply(lambda x: x.get('Note', ''))

    grouped_exon_gff_df = exon_gff_df.groupby(unique_cols).aggregate(list)

    # Get existing exon mappings using SQLAlchemy 2.0 style
    stmt = select(base.Exon)
    vdb_exon_idx_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id 
        for e in session.execute(stmt).scalars().all()
    }

    exons = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    for exon_idx, row in grouped_exon_gff_df.iterrows():
        if exon_idx in vdb_exon_idx_map.keys():
            exon_id = vdb_exon_idx_map[exon_idx]
            
            update_entries.append({
                "id": exon_id, 
                "refseq_id": row.ID[0],
                "attrs": {
                    "product": row['product'][0] if 'product' in row else '',
                    "parent": row.Parent[0]
                }
            })
            
        elif exon_idx not in existing_entries:
            existing_entries.append(exon_idx)
        
            exons.append(base.Exon(
                start=exon_idx[0],
                end=exon_idx[1],
                strand=exon_idx[2],
                assembly_id=exon_idx[3],
                ensembl_id='',
                refseq_id=row.ID[0],
                chess_id='',
                velia_id='',
                attrs={
                    "inference": row.inference[0] if 'inference' in row else '',
                    "model_evidence": row.model_evidence[0] if 'model_evidence' in row else '',
                    "notes": row.Note[0]
                }
            ))
        
        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            if row.Name[i]:
                syn_dict[row.source[i]].add(row.Name[i])
            if ID:
                syn_dict[row.source[i]].add(ID)

        synonym_dict[exon_idx] = syn_dict

    # Bulk update existing exons
    if update_entries:
        session.bulk_update_mappings(base.Exon, update_entries)
        session.commit()

    # Bulk insert new exons
    if exons:
        session.bulk_save_objects(exons)
        session.commit()

    logging.info(f'Updated {len(update_entries)} exons with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(exons)} RefSeq exons without an exact coordinate match to GENCODE')

    # Get exon IDs for xrefs using SQLAlchemy 2.0 style
    stmt = select(base.Exon)
    vdb_exon_idx_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id 
        for e in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs using SQLAlchemy 2.0 style
    stmt = select(base.Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }

    # Prepare xrefs for bulk insert
    exon_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                exon_xrefs.append(base.SequenceRegionXref(
                    sequence_region_id=vdb_exon_idx_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['RefSeq'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    # Bulk insert xrefs
    if exon_xrefs:
        session.bulk_save_objects(exon_xrefs)
        session.commit()

    logging.info(f'Added {len(exon_xrefs)} RefSeq exon synonyms')


def load_refseq_cds(session: Session, cds_gff_df: pd.DataFrame) -> None:
    """Load CDS from RefSeq GFF file.

    Args:
        session: SQLAlchemy session object
        cds_gff_df: RefSeq GFF dataframe filtered for CDS records
    """
    # Group by unique columns to handle redundant entries
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    cds_gff_df.sort_values(by=unique_cols, inplace=True)

    # Process attributes
    cds_gff_df['attrs'] = cds_gff_df['attributes'].apply(parse_attributes)
    cds_gff_df['ID'] = cds_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    cds_gff_df['Name'] = cds_gff_df['attrs'].apply(lambda x: x.get('Name', ''))
    cds_gff_df['protein_id'] = cds_gff_df['attrs'].apply(lambda x: x.get('protein_id', ''))

    grouped_cds_gff_df = cds_gff_df.groupby(unique_cols).aggregate(list)

    # Get existing CDS mappings using SQLAlchemy 2.0 style
    stmt = select(base.Cds)
    vdb_cds_idx_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id 
        for c in session.execute(stmt).scalars().all()
    }

    cds = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    for cds_idx, row in grouped_cds_gff_df.iterrows():
        if cds_idx in vdb_cds_idx_map.keys():
            cds_id = vdb_cds_idx_map[cds_idx]
            
            update_entries.append({
                "id": cds_id, 
                "refseq_id": row.ID[0], 
                "refseq_protein_id": row.protein_id[0]
            })
            
        elif cds_idx not in existing_entries:
            existing_entries.append(cds_idx)
        
            cds.append(base.Cds(
                start=cds_idx[0],
                end=cds_idx[1],
                strand=cds_idx[2],
                assembly_id=cds_idx[3],
                ensembl_id='',
                refseq_id=row.ID[0],
                chess_id='',
                velia_id='',
                refseq_protein_id=row.protein_id[0]
            ))

        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            if row.Name[i]:
                syn_dict[row.source[i]].add(row.Name[i])
            if ID:
                syn_dict[row.source[i]].add(ID)

        synonym_dict[cds_idx] = syn_dict

    # Bulk update existing CDS
    if update_entries:
        session.bulk_update_mappings(base.Cds, update_entries)
        session.commit()
    
    # Bulk insert new CDS
    if cds:
        session.bulk_save_objects(cds)
        session.commit()

    logging.info(f'Updated {len(update_entries)} CDS with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(cds)} RefSeq CDS without an exact coordinate match to GENCODE')

    # Get CDS IDs for xrefs using SQLAlchemy 2.0 style
    stmt = select(base.Cds)
    vdb_cds_idx_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id 
        for c in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs using SQLAlchemy 2.0 style
    stmt = select(base.Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }

    # Prepare xrefs for bulk insert
    cds_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                cds_xrefs.append(base.SequenceRegionXref(
                    sequence_region_id=vdb_cds_idx_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['RefSeq'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    # Bulk insert xrefs
    if cds_xrefs:
        session.bulk_save_objects(cds_xrefs)
        session.commit()

    logging.info(f'Added {len(cds_xrefs)} RefSeq CDS synonyms')


def load_refseq_lncRNAs(session: Session, lnc_gff_df: pd.DataFrame) -> None:
    """Load lncRNA transcripts from RefSeq GFF file.

    This function is almost identical to load transcripts. In the GENCODE loading
    there is no distinction for transcripts/lncRNAs. However, in the RefSeq GFF
    there are a few different fields and clearly a different annotation process that
    make lncRNAs look more like genes. I still load them here as transcripts but
    it should be possible to merge this code with transcript loading in the future.

    Args:
        session: SQLAlchemy session object
        lnc_gff_df: RefSeq GFF dataframe filtered for lncRNA records
    """
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    lnc_gff_df.sort_values(by=unique_cols, inplace=True)

    # Process attributes
    lnc_gff_df['attrs'] = lnc_gff_df['attributes'].apply(parse_attributes)
    lnc_gff_df['ID'] = lnc_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    lnc_gff_df['Name'] = lnc_gff_df['attrs'].apply(lambda x: x.get('Name', ''))
    lnc_gff_df['transcript_id'] = lnc_gff_df['attrs'].apply(lambda x: x.get('transcript_id', ''))
    lnc_gff_df['entrez_gene_id'] = lnc_gff_df['attrs'].apply(lambda x: x.get('entrez_gene_id', ''))
    lnc_gff_df['gbkey'] = lnc_gff_df['attrs'].apply(lambda x: x.get('gbkey', ''))
    lnc_gff_df['experiment'] = lnc_gff_df['attrs'].apply(lambda x: x.get('experiment', ''))
    lnc_gff_df['model_evidence'] = lnc_gff_df['attrs'].apply(lambda x: x.get('model_evidence', ''))
    lnc_gff_df['product'] = lnc_gff_df['attrs'].apply(lambda x: x.get('product', ''))

    grouped_lnc_gff_df = lnc_gff_df.groupby(unique_cols).aggregate(list)

    # Get gene mappings using SQLAlchemy 2.0 style
    stmt = select(SequenceRegionXref).join(Gene, Gene.id == SequenceRegionXref.sequence_region_id)
    synonym_vdb_gene_map = {gs.xref: gs.sequence_region_id for gs in session.execute(stmt).scalars().all()}

    # Get transcript mappings using SQLAlchemy 2.0 style
    stmt = select(Transcript)
    vdb_tx_idx_map = {(t.start, t.end, t.strand, t.assembly_id): t.id for t in session.execute(stmt).scalars().all()}

    transcripts = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    # Get unmapped gene using SQLAlchemy 2.0 style
    stmt = select(Gene).filter_by(hgnc_id='unmapped')
    unmapped_gene = session.execute(stmt).scalar_one()

    for tx_idx, row in grouped_lnc_gff_df.iterrows():
        try:
            gene_id = synonym_vdb_gene_map[row.entrez_gene_id[0]]
        except:
            gene_id = unmapped_gene.id

        if tx_idx in vdb_tx_idx_map.keys():
            tx_id = vdb_tx_idx_map[tx_idx]

            update_entries.append({"id": tx_id, "refseq_id": row.transcript_id[0],
                                "attrs": {"experiment": row.experiment[0],
                                          "model_evidence": row.model_evidence[0],
                                          "product": row['product'][0]}})

        elif tx_idx not in existing_entries:
            existing_entries.append(tx_idx)

            transcripts.append(Transcript(tx_idx[0], tx_idx[1], tx_idx[2], tx_idx[3],
                                        gene_id, '', row.Name[0], '', '', 
                                        '', '', '', row.gbkey[0], 
                                        attrs={"experiment": row.experiment[0],
                                               "model_evidence": row.model_evidence[0],
                                               "product": row['product'][0]}))

        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.Name[i])
            syn_dict[row.source[i]].add(ID)
            syn_dict[row.source[i]].add(row.transcript_id[i])

        synonym_dict[tx_idx] = syn_dict

    session.bulk_update_mappings(Transcript, update_entries)
    session.commit()

    session.add_all(transcripts)
    session.commit()

    logging.info(f'Updated {len(update_entries)} lncRNAs with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(transcripts)} RefSeq lncRNAs without an exact coordinate match to GENCODE')

    vdb_tx_idx_map = {(t.start, t.end, t.strand, t.assembly_id): t.id for t in session.query(Transcript).all()}

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    transcript_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '': continue
                transcript_xrefs.append(SequenceRegionXref(
                    sequence_region_id=vdb_tx_idx_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['RefSeq'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    session.add_all(transcript_xrefs)
    session.commit()


def load_refseq_transcript_exons(session: Session, transcript_exon_map: Dict[str, List[Tuple[int, int]]]) -> None:
    """Load transcript-exon relationships from RefSeq data.

    Args:
        session: SQLAlchemy session object
        transcript_exon_map: Mapping of transcript IDs to lists of (exon_number, exon_id) tuples
    """
    # Get RefSeq dataset ID using SQLAlchemy 2.0 style
    stmt = select(Dataset.id).filter(Dataset.name == 'RefSeq')
    refseq_dataset_id = session.execute(stmt).scalar_one()

    # Get transcript mappings using SQLAlchemy 2.0 style
    stmt = select(Transcript)
    refseq_vdb_transcript_map = {
        t.refseq_id: t.id 
        for t in session.execute(stmt).scalars().all()
    }

    # Get additional transcript mappings from xrefs
    stmt = select(TranscriptXref).filter(
        and_(
            TranscriptXref.transcript_dataset_id == refseq_dataset_id,
            TranscriptXref.xref_dataset_id == refseq_dataset_id
        )
    )
    refseq_vdb_transcript_map.update({
        tx.xref: tx.transcript_id 
        for tx in session.execute(stmt).scalars().all()
    })
    
    transcript_exons = []

    # Get existing relationships using SQLAlchemy 2.0 style
    stmt = select(TranscriptExon)
    existing_entries = {
        (te.transcript_id, te.exon_id, te.exon_number) 
        for te in session.execute(stmt).scalars().all()
    }
    entries = set()

    for transcript_id, vals in transcript_exon_map.items():
        for (exon_number, exon_id) in vals:
            entries.add((refseq_vdb_transcript_map[transcript_id], exon_id, exon_number))

    entries.difference_update(existing_entries)
            
    for entry in entries:
        transcript_exons.append(TranscriptExon(
            transcript_id=entry[0],
            exon_id=entry[1],
            exon_number=entry[2],
            attrs={'mapping_source': 'RefSeq'}
        ))
        
    session.bulk_save_objects(transcript_exons)
    session.commit()

    logging.info(f'Added {len(transcript_exons)} transcript <-> exon mappings from RefSeq')


def load_refseq_exon_cds():
    """
    TODO: Additional code needs to be written to try and map
    RefSeq exons to RefSeq CDS.  For some reason this mapping
    is not provided in the GFF and so a position based lookup
    would be needed to enable loading this for RefSeq
    """
    return


def load_chess_exons(session: Session, exon_gff_df: pd.DataFrame) -> None:
    """Load exons from CHESS GFF file.

    Args:
        session: SQLAlchemy session object
        exon_gff_df: CHESS GFF dataframe filtered for exon records
    """
    # Group by unique columns to handle redundant entries
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    exon_gff_df.sort_values(by=unique_cols, inplace=True)

    # Process attributes
    exon_gff_df['attrs'] = exon_gff_df['attributes'].apply(parse_attributes)
    exon_gff_df['ID'] = exon_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    exon_gff_df['Name'] = exon_gff_df['attrs'].apply(lambda x: x.get('Name', ''))
    exon_gff_df['Parent'] = exon_gff_df['attrs'].apply(lambda x: x.get('Parent', ''))
    exon_gff_df['gene_id'] = exon_gff_df['attrs'].apply(lambda x: x.get('gene_id', ''))
    exon_gff_df['transcript_id'] = exon_gff_df['attrs'].apply(lambda x: x.get('transcript_id', ''))

    grouped_exon_gff_df = exon_gff_df.groupby(unique_cols).aggregate(list)

    # Get existing exon mappings using SQLAlchemy 2.0 style
    stmt = select(base.Exon)
    vdb_exon_idx_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id 
        for e in session.execute(stmt).scalars().all()
    }

    exons = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    for exon_idx, row in grouped_exon_gff_df.iterrows():
        if exon_idx in vdb_exon_idx_map.keys():
            exon_id = vdb_exon_idx_map[exon_idx]
            
            update_entries.append({
                "id": exon_id, 
                "chess_id": row.ID[0],
                "attrs": {
                    "parent": row.Parent[0],
                    "gene_id": row.gene_id[0],
                    "transcript_id": row.transcript_id[0]
                }
            })
            
        elif exon_idx not in existing_entries:
            existing_entries.append(exon_idx)
        
            exons.append(base.Exon(
                start=exon_idx[0],
                end=exon_idx[1],
                strand=exon_idx[2],
                assembly_id=exon_idx[3],
                ensembl_id='',
                refseq_id='',
                chess_id=row.ID[0],
                velia_id='',
                attrs={
                    "parent": row.Parent[0],
                    "gene_id": row.gene_id[0],
                    "transcript_id": row.transcript_id[0]
                }
            ))
        
        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            if row.Name[i]:
                syn_dict[row.source[i]].add(row.Name[i])
            if ID:
                syn_dict[row.source[i]].add(ID)

        synonym_dict[exon_idx] = syn_dict

    # Bulk update existing exons
    if update_entries:
        session.bulk_update_mappings(base.Exon, update_entries)
        session.commit()

    # Bulk insert new exons
    if exons:
        session.bulk_save_objects(exons)
        session.commit()

    logging.info(f'Updated {len(update_entries)} exons with CHESS info that had exact coordinate matches')
    logging.info(f'Added {len(exons)} CHESS exons without an exact coordinate match')

    # Get exon IDs for xrefs using SQLAlchemy 2.0 style
    stmt = select(base.Exon)
    vdb_exon_idx_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id 
        for e in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs using SQLAlchemy 2.0 style
    stmt = select(base.Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }

    # Prepare xrefs for bulk insert
    exon_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                exon_xrefs.append(base.SequenceRegionXref(
                    sequence_region_id=vdb_exon_idx_map[idx],
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['CHESS'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    # Bulk insert xrefs
    if exon_xrefs:
        session.bulk_save_objects(exon_xrefs)
        session.commit()

    logging.info(f'Added {len(exon_xrefs)} CHESS exon synonyms')


def load_chess_transcripts(session: Session, tx_gff_df: pd.DataFrame, exon_gff_df: pd.DataFrame) -> Dict[str, List[Tuple[int, int]]]:
    """Load transcripts from CHESS GFF file.

    Args:
        session: SQLAlchemy session object
        tx_gff_df: CHESS GFF dataframe filtered for transcript records
        exon_gff_df: CHESS GFF dataframe filtered for exon records

    Returns:
        Dict[str, List[Tuple[int, int]]]: Mapping of transcript IDs to lists of (exon_number, exon_id) tuples
    """
    exon_gff_df.sort_values(by='start', inplace=True)

    # Process attributes for both dataframes
    exon_gff_df['attrs'] = exon_gff_df['attributes'].apply(parse_attributes)
    exon_gff_df['ID'] = exon_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    exon_gff_df['Parent'] = exon_gff_df['attrs'].apply(lambda x: x.get('Parent', ''))
    exon_gff_df['gene_id'] = exon_gff_df['attrs'].apply(lambda x: x.get('gene_id', ''))
    exon_gff_df['transcript_id'] = exon_gff_df['attrs'].apply(lambda x: x.get('transcript_id', ''))

    tx_gff_df['attrs'] = tx_gff_df['attributes'].apply(parse_attributes)
    tx_gff_df['ID'] = tx_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    tx_gff_df['geneID'] = tx_gff_df['attrs'].apply(lambda x: x.get('geneID', ''))
    tx_gff_df['gene_name'] = tx_gff_df['attrs'].apply(lambda x: x.get('gene_name', ''))
    tx_gff_df['gene_type'] = tx_gff_df['attrs'].apply(lambda x: x.get('gene_type', ''))
    tx_gff_df['db_xref'] = tx_gff_df['attrs'].apply(lambda x: x.get('db_xref', ''))
    tx_gff_df['num_samples'] = tx_gff_df['attrs'].apply(lambda x: x.get('num_samples', ''))
    tx_gff_df['max_tpm'] = tx_gff_df['attrs'].apply(lambda x: x.get('max_tpm', ''))

    exon_grouped_df = exon_gff_df.groupby('Parent').aggregate(list)

    exon_grouped_df['exon_number'] = exon_grouped_df.apply(
        lambda x: range(1, len(x.start)+1) if x.strand[0] == '+' 
        else [y for y in range(1, len(x.start)+1)][::-1], 
        axis=1
    )
    
    transcript_exon_gff_df = exon_grouped_df.merge(
        tx_gff_df, 
        left_index=True, 
        right_on='ID', 
        suffixes=('_ex', '_tx'), 
        how='left'
    )
    transcript_exon_gff_df.set_index(['ID', 'start_tx', 'end_tx', 'strand_tx', 'assembly_id_tx'], inplace=True)

    # Get gene mappings using SQLAlchemy 2.0 style
    stmt = select(SequenceRegionXref).join(Gene, Gene.id == SequenceRegionXref.sequence_region_id)
    synonym_vdb_gene_map = {
        gs.xref: gs.sequence_region_id 
        for gs in session.execute(stmt).scalars().all()
    }

    # Get transcript and exon mappings using SQLAlchemy 2.0 style
    stmt = select(Transcript)
    vdb_transcript_idx_map = {
        t.transcript_idx: t.id 
        for t in session.execute(stmt).scalars().all()
    }

    stmt = select(Exon)
    vdb_exon_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id 
        for e in session.execute(stmt).scalars().all()
    }

    transcripts = []
    transcript_exon_map = {}
    synonym_dict = {}

    existing_entries = []
    update_entries = []
    gene_update_entries = []

    # Get unmapped gene using SQLAlchemy 2.0 style
    stmt = select(Gene).filter_by(hgnc_id='unmapped')
    unmapped_gene = session.execute(stmt).scalar_one()

    for (transcript_id, start, end, strand, assembly_id), row in transcript_exon_gff_df.iterrows():
        try:
            gene_id = synonym_vdb_gene_map[row.gene_name]
            gene_update_entries.append({"id": gene_id, "chess_id": row.geneID})
        except:
            gene_id = unmapped_gene.id
        
        block_sizes = []
        chrom_starts = []
        exon_ids = []

        for i in range(len(row.ID_ex)):
            if pd.isna(row.ID_ex[i]):
                continue

            chrom_starts.append(str(row.start_ex[i]))
            block_sizes.append(str(row.end_ex[i] - row.start_ex[i]))
            exon_ids.append((row.exon_number[i], vdb_exon_map[(row.start_ex[i], row.end_ex[i], row.strand_ex[i], assembly_id)]))

        block_sizes = ';'.join(block_sizes)
        chrom_starts = ';'.join(chrom_starts)
        transcript_idx_str = f'{start}_{end}_{strand}_{assembly_id}_{block_sizes}_{chrom_starts}'
        transcript_idx = hashlib.sha3_512(transcript_idx_str.encode('utf-8')).hexdigest()

        if transcript_idx in vdb_transcript_idx_map.keys():
            transcript_vdb_id = vdb_transcript_idx_map[transcript_idx]

            update_entries.append({
                "id": transcript_vdb_id, 
                "chess_id": transcript_id, 
                "attrs": {
                    "gene_type": row.gene_type,
                    "db_xref": row.db_xref,
                    "num_samples": row.num_samples,
                    "max_tpm": row.max_tpm
                }
            })

        elif transcript_idx not in existing_entries:
            existing_entries.append(transcript_idx)
            
            transcripts.append(Transcript(
                start=start,
                end=end,
                strand=strand,
                assembly_id=assembly_id,
                block_sizes=block_sizes,
                chrom_starts=chrom_starts,
                transcript_idx=transcript_idx,
                transcript_idx_str=transcript_idx_str,
                gene_id=gene_id,
                ensembl_id='',
                refseq_id='',
                chess_id=transcript_id,
                velia_id='',
                support_level='',
                transcript_type=row.gene_type,
                attrs={
                    "db_xref": row.db_xref,
                    "num_samples": row.num_samples,
                    "max_tpm": row.max_tpm
                }
            ))

        syn_dict = {source: set() for source in row.source_ex}

        for i, ID in enumerate(row.start_ex):
            syn_dict[row.source_ex[i]].add(row.db_xref)
            syn_dict[row.source_ex[i]].add(transcript_id)

        synonym_dict[transcript_idx] = syn_dict

        transcript_exon_map[transcript_idx] = []
        for (exon_num, exon_id) in exon_ids:
            transcript_exon_map[transcript_idx].append((exon_num, exon_id))

    # Bulk update existing records
    if gene_update_entries:
        session.bulk_update_mappings(Gene, gene_update_entries)
        session.commit()

    if update_entries:
        session.bulk_update_mappings(Transcript, update_entries)
        session.commit()

    # Bulk insert new transcripts
    if transcripts:
        session.bulk_save_objects(transcripts)
        session.commit()

    logging.info(f'Updated {len(gene_update_entries)} genes with CHESS info that had a synonym match to GENCODE/RefSeq')
    logging.info(f'Updated {len(update_entries)} transcripts with CHESS info that had exact coordinate matches to GENCODE/RefSeq')
    logging.info(f'Added {len(transcripts)} CHESS transcripts without an exact coordinate match to GENCODE/RefSeq')

    # Get transcript IDs for xrefs using SQLAlchemy 2.0 style
    stmt = select(Transcript)
    vdb_transcript_idx_map = {
        t.transcript_idx: t.id 
        for t in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs using SQLAlchemy 2.0 style
    stmt = select(Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }

    transcript_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                transcript_xrefs.append(TranscriptXref(
                    transcript_id=vdb_transcript_idx_map[idx],
                    xref=synonym,
                    type='synonym',
                    transcript_dataset_id=dataset_ids['CHESS'],
                    xref_dataset_id=dataset_ids['CHESS']
                ))

    # Bulk insert xrefs
    if transcript_xrefs:
        session.bulk_save_objects(transcript_xrefs)
        session.commit()

    return transcript_exon_map


def load_chess_transcript_exons(session: Session, transcript_exon_map: Dict[str, List[Tuple[int, int]]]) -> None:
    """Load transcript-exon relationships from CHESS data.

    Args:
        session: SQLAlchemy session object
        transcript_exon_map: Mapping of transcript IDs to lists of (exon_number, exon_id) tuples
    """
    # Get transcript mappings using SQLAlchemy 2.0 style
    stmt = select(Transcript)
    vdb_transcript_map = {
        t.transcript_idx: t.id 
        for t in session.execute(stmt).scalars().all()
    }
    
    transcript_exons = []

    stmt = select(TranscriptExon)
    existing_entries = {
        (te.transcript_id, te.exon_id, te.exon_number) 
        for te in session.execute(stmt).scalars().all()
    }
    entries = set()

    for transcript_idx, vals in transcript_exon_map.items():
        for (exon_number, exon_id) in vals:
            entries.add((vdb_transcript_map[transcript_idx], exon_id, exon_number))

    entries.difference_update(existing_entries)
            
    for entry in entries:
        transcript_exons.append(TranscriptExon(
            transcript_id=entry[0],
            exon_id=entry[1],
            exon_number=entry[2],
            attrs={'mapping_source': 'CHESS'}
        ))
        
    if transcript_exons:
        session.bulk_save_objects(transcript_exons)
        session.commit()

    logging.info(f'Added {len(transcript_exons)} transcript <-> exon mappings from CHESS')


def load_chess_cds(session: Session, cds_gff_df: pd.DataFrame) -> None:
    """Load CDS data from CHESS database.

    This function processes and loads coding sequence (CDS) annotations from CHESS,
    creating new CDS entries and their cross-references.

    Args:
        session: SQLAlchemy session object
        cds_gff_df: CHESS GFF dataframe filtered for CDS regions
    """
    # Group by unique columns to handle redundant entries
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    cds_gff_df.sort_values(by=unique_cols, inplace=True)

    # Process attributes
    cds_gff_df['attrs'] = cds_gff_df['attributes'].apply(parse_attributes)
    cds_gff_df['ID'] = cds_gff_df['attrs'].apply(lambda x: x.get('ID', ''))
    cds_gff_df['Parent'] = cds_gff_df['attrs'].apply(lambda x: x.get('Parent', ''))
    cds_gff_df['db_xref'] = cds_gff_df['attrs'].apply(lambda x: x.get('db_xref', ''))

    grouped_cds_gff_df = cds_gff_df.groupby(unique_cols).aggregate(list)

    # Get existing CDS mappings using SQLAlchemy 2.0 style
    stmt = select(Cds)
    vdb_cds_idx_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id 
        for c in session.execute(stmt).scalars().all()
    }

    cds = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    for cds_idx, row in grouped_cds_gff_df.iterrows():
        chess_id = f'cds-{row.Parent[0]}'

        if cds_idx in vdb_cds_idx_map.keys():
            cds_id = vdb_cds_idx_map[cds_idx]
            update_entries.append({
                "id": cds_id, 
                "chess_id": chess_id
            })

        elif cds_idx not in existing_entries:
            existing_entries.append(cds_idx)
            
            cds.append(Cds(
                start=cds_idx[0],
                end=cds_idx[1],
                strand=cds_idx[2],
                assembly_id=cds_idx[3],
                ensembl_id='',
                refseq_id='',
                chess_id=chess_id,
                velia_id='',    
                attrs={
                    "db_xref": row.db_xref
                }
            ))

        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(f'cds-{row.Parent[i]}')
            syn_dict[row.source[i]].add(f'cds-{row.db_xref[i]}')

        synonym_dict[cds_idx] = syn_dict

    # Bulk update existing CDS
    if update_entries:
        session.bulk_update_mappings(Cds, update_entries)
        session.commit()
    
    # Bulk insert new CDS
    if cds:
        session.bulk_save_objects(cds)
        session.commit()

    logging.info(f'Updated {len(update_entries)} CDS with CHESS info that had exact coordinate matches to GENCODE/RefSeq')
    logging.info(f'Added {len(cds)} CHESS CDS without an exact coordinate match to GENCODE/RefSeq')

    # Get CDS IDs for xrefs using SQLAlchemy 2.0 style
    stmt = select(Cds)
    vdb_cds_idx_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id 
        for c in session.execute(stmt).scalars().all()
    }

    # Get dataset IDs using SQLAlchemy 2.0 style
    stmt = select(Dataset)
    dataset_ids = {
        d.name: d.id 
        for d in session.execute(stmt).scalars().all()
    }

    stmt = select(SequenceRegionXref)
    existing_xrefs = {
        (x.sequence_region_id, x.xref, x.type, x.sequence_region_dataset_id, x.xref_dataset_id)
        for x in session.execute(stmt).scalars().all()
    }

    # Prepare xrefs for bulk insert
    cds_xrefs = []
    for idx, synonym_entries in synonym_dict.items():
        cds_id = vdb_cds_idx_map[idx]
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                    
                xref_tuple = (
                    cds_id,
                    synonym,
                    'synonym',
                    dataset_ids['CHESS'],
                    dataset_ids[dataset_name]
                )
                
                # Skip if this combination already exists
                if xref_tuple in existing_xrefs:
                    continue
                    
                cds_xrefs.append(SequenceRegionXref(
                    sequence_region_id=cds_id,
                    xref=synonym,
                    type='synonym',
                    sequence_region_dataset_id=dataset_ids['CHESS'],
                    xref_dataset_id=dataset_ids[dataset_name]
                ))

    # Bulk insert xrefs
    if cds_xrefs:
        session.bulk_save_objects(cds_xrefs)
        session.commit()

    logging.info(f'Added {len(cds_xrefs)} CHESS CDS synonyms')


def load_gencode_lncRNA_genes(session, gene_gff_df):
    """Load long non-coding RNA genes from GENCODE GFF data.

    Args:
        session: SQLAlchemy session object
        gene_gff_df (pd.DataFrame): GENCODE GFF dataframe filtered for lncRNA genes

    Note:
        Creates both gene entries and cross-references in the database.
        Handles duplicates by updating existing records.
    """
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    gene_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_gene_gff_df = gene_gff_df.groupby(unique_cols).aggregate(list)

    vdb_gene_idx_map = {
        (g.start, g.end, g.strand, g.assembly_id): g.id 
        for g in session.query(Gene).all()
    }

    update_entries = []
    new_genes = []
    synonym_dict = {}

    for gene_idx, row in grouped_gene_gff_df.iterrows():
        if gene_idx in vdb_gene_idx_map:
            gene_id = vdb_gene_idx_map[gene_idx]
            update_entries.append({
                "id": gene_id,
                "ensembl_id": row.ID[0],
                "attrs": {
                    "gene_type": row.gene_type[0],
                    "tag": row.tag[0]
                }
            })
        elif gene_idx not in synonym_dict:
            new_genes.append(Gene(
                gene_idx[0], gene_idx[1], gene_idx[2], gene_idx[3],
                row.hgnc_id[0], row.gene_name[0], row.ID[0], 
                '', '', '', row.gene_type[0], '',
                attrs={
                    "gene_type": row.gene_type[0], 
                    "tag": row.tag[0]
                }
            ))

        syn_dict = {source: set() for source in row.source}
        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.Name[i])
            syn_dict[row.source[i]].add(ID)
            syn_dict[row.source[i]].add(row.gene_id[i])

        syn_dict['hgnc_id'] = set(row.hgnc_id)
        synonym_dict[gene_idx] = syn_dict

    session.bulk_update_mappings(Gene, update_entries)
    session.commit()

    session.add_all(new_genes)
    session.commit()

    logging.info(
        f'Updated {len(update_entries)} GENCODE lncRNA genes that had exact '
        'coordinate matches to Refseq'
    )
    logging.info(f'Added {len(new_genes)} GENCODE lncRNA genes')

    gene_xrefs = []
    vdb_gene_idx_map = {
        (g.start, g.end, g.strand, g.assembly_id): g.id
        for g in session.query(Gene).all()
    }
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if not synonym:
                    continue
                gene_xrefs.append(SequenceRegionXref(
                    vdb_gene_idx_map[idx],
                    synonym,
                    'synonym',
                    dataset_ids['ENSEMBL'],
                    dataset_ids[dataset_name]
                ))

    session.add_all(gene_xrefs)
    session.commit()

    logging.info(f'Added {len(gene_xrefs)} GENCODE lncRNA synonyms')


def load_bigprot_orfs(session: Session, bigprot_dir: Path) -> None:
    """Load BigProt ORF entries from a CSV file.

    Args:
        session: SQLAlchemy session object
        bigprot_dir: Directory containing BigProt files
    """
    chunk_size = 600000
    bigprot_csv = bigprot_dir.joinpath(
        'orfset_BigProt_minlen_15_maxlen_999999999_orfs.csv.gz'
    )

    reader = pd.read_csv(bigprot_csv, chunksize=chunk_size)

    for chunk in reader:
        orfs = []
        chunk.columns = [col[4:] for col in chunk.columns]

        for _, row in chunk.iterrows():
            orfs.append(Orf(
                start=row.start,
                end=row.end,
                strand=row.strand,
                assembly_id=row.assembly_id,
                block_sizes=row.block_sizes,
                chrom_starts=row.chrom_starts,
                phases=row.phases,
                reading_frames=row.exon_frames,
                orf_idx=row.orf_idx,
                orf_idx_str=row.orf_idx_str,
                secondary_orf_id=row.secondary_orf_id,
                aa_seq=row.aa_seq,
                id=row.id
            ))
        
        session.bulk_save_objects(orfs)
        session.commit()


def load_bigprot_transcript_orfs(session: Session, bigprot_dir: Path) -> None:
    """Load BigProt transcript entries from a CSV file.

    Args:
        session: SQLAlchemy session object
        bigprot_dir: Directory containing BigProt files
    """
    chunk_size = 600000
    bigprot_csv = bigprot_dir.joinpath(
        f'orfset_BigProt_minlen_15_maxlen_999999999_transcript_orfs.csv.gz'
    )

    reader = pd.read_csv(bigprot_csv, chunksize=chunk_size)

    for chunk in reader:
        tx_orfs = []

        for _, row in chunk.iterrows():
            tx_orfs.append(TranscriptOrf(
                transcript_id=row['transcript_orf.transcript_id'],
                orf_id=row['transcript_orf.orf_id'],
                evidence_tag=row['transcript_orf.evidence_tag']
            ))
        
        session.bulk_save_objects(tx_orfs)
        session.commit()


def load_bigprot_cds_orf(session: Session, bigprot_dir: Path) -> None:
    """Load BigProt CDS ORF entries from a CSV file.

    Args:
        session: SQLAlchemy session object
        bigprot_dir: Directory containing BigProt files
    """
    chunk_size = 500000
    bigprot_csv = bigprot_dir.joinpath(
        f'orfset_BigProt_minlen_15_maxlen_999999999_cds_orfs.csv.gz'
    )

    reader = pd.read_csv(bigprot_csv, chunksize=chunk_size)

    for chunk in reader:
        cds_orfs = []

        for _, row in chunk.iterrows():
            cds_orfs.append(CdsOrf(
                cds_id=row['cds_orf.cds_id'],
                orf_id=row['cds_orf.orf_id'],
                cds_number=row['cds_orf.cds_number'],
                phase=row['cds_orf.phase'],
                reading_frame=row['cds_orf.reading_frame']
            ))
        
        session.bulk_save_objects(cds_orfs)
        session.commit()


def parse_attributes(attr_str: str) -> Dict[str, Any]:
    """Parse GFF/GTF attribute string into a dictionary.
    
    Args:
        attr_str: String containing semicolon-separated key-value pairs
        
    Returns:
        Dictionary of attribute key-value pairs
        
    Example:
        >>> parse_attributes('ID=gene1;Name=BRCA1;gene_type=protein_coding')
        {'ID': 'gene1', 'Name': 'BRCA1', 'gene_type': 'protein_coding'}
    """
    if not attr_str or pd.isna(attr_str):
        return {}
        
    attrs = {}
    for attr in attr_str.split(';'):
        if not attr.strip():
            continue
            
        try:
            key, value = attr.strip().split('=', 1)
            attrs[key.strip()] = value.strip(' "\'')
        except ValueError:
            # Skip malformed attributes
            continue
            
    return attrs
