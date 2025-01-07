"""Module to update core annotation data"""

from Bio.Seq import Seq
from orfdb.base import *
from orfdb import util as orf_utils

import ast
import gzip
import hashlib
import logging
import pandas as pd
from os.path import basename
from warnings import warn


def update_ensembl_entrez_gene_mapping(gencode_dir, version, session):
    """Update Ensembl to Entrez gene ID mappings.

    Args:
        gencode_dir (Path): Directory containing GENCODE files
        version (str): GENCODE version
        session: SQLAlchemy session object
    """
    entrez_gencode_df = pd.read_csv(
        gencode_dir.joinpath(
            version,
            f'gencode.{version}.metadata.EntrezGene.gz'
        ),
        sep='\t',
        names=['ensembl_id', 'entrez_id']
    )

    gene_tx_map = session.query(Gene.ensembl_id, Transcript.ensembl_id)\
        .join(Transcript, Transcript.gene_id == Gene.id)\
        .filter(Transcript.ensembl_id != '').all()

    ensembl_map_df = pd.DataFrame(
        gene_tx_map,
        columns=['ensembl_gene_id', 'ensembl_transcript_id']
    )

    merged_df = entrez_gencode_df.merge(
        ensembl_map_df,
        left_on='ensembl_id',
        right_on='ensembl_transcript_id'
    )

    entrez_gencode_map_df = merged_df[['entrez_id', 'ensembl_gene_id']]\
        .drop_duplicates()

    ensembl_orf_gene_map = {
        g.ensembl_id: (g.id, g.refseq_id)
        for g in session.query(Gene).all()
    }
    
    update_entries = []
    synonym_dict = {}

    for i, row in entrez_gencode_map_df.iterrows():
        gene_id = ensembl_orf_gene_map[row.ensembl_gene_id]
        
        if gene_id[1] != '' and gene_id[1] != str(row.entrez_id):
            synonym_dict[gene_id[0]] = {'RefSeq': [row.entrez_id]}
        
        update_entries.append({
            "id": gene_id[0],
            "refseq_id": row.entrez_id
        })

    session.bulk_update_mappings(Gene, update_entries)
    session.commit()

    logging.info(
        f'Updated {len(update_entries)} genes with EntrezID that are mapped '
        'to GENCODE gene'
    )

    gene_xrefs = []
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for gene_id, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '':
                    continue
                gene_xrefs.append(
                    SequenceRegionXref(
                        gene_id,
                        synonym,
                        'synonym',
                        dataset_ids['RefSeq'],
                        dataset_ids[dataset_name]
                    )
                )

    session.add_all(gene_xrefs)
    session.commit()

    logging.info(f'Added {len(gene_xrefs)} RefSeq gene synonyms')


def update_benchling_protein_xrefs(session):
    """Update Benchling protein cross-references.

    Args:
        session: SQLAlchemy session object
    """
    pass


def update_chess_transcript_ids(chess_dir, session):
    """Update CHESS transcript IDs.

    Args:
        chess_dir (Path): Directory containing CHESS files
        session: SQLAlchemy session object
    """
    chess_expanded_gff = chess_dir.joinpath('chess3.0.expanded.gff')
    
    chess_df = pd.read_csv(chess_expanded_gff, sep='\t', low_memory=False)
    chess_df.fillna('', inplace=True)
    
    tx_gff_df = chess_df[chess_df['type'] == 'transcript'].copy()
    exon_gff_df = chess_df[chess_df['type'] == 'exon'].copy()
    
    exon_gff_df.sort_values(by='start', inplace=True)
    exon_grouped_df = exon_gff_df.groupby('Parent').aggregate(list)
    exon_grouped_df['exon_number'] = exon_grouped_df.apply(
        lambda x: range(1, len(x.start)+1)
        if x.strand[0] == '+'
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
    transcript_exon_gff_df.set_index(
        ['ID', 'start_tx', 'end_tx', 'strand_tx', 'assembly_id_tx'],
        inplace=True
    )

    orf_exon_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id
        for e in session.query(Exon).all()
    }
    orf_transcript_idx_map = {
        t.transcript_idx: t.id
        for t in session.query(Transcript).all()
    }

    update_entries = []

    for (transcript_id, start, end, strand, assembly_id), row in \
            transcript_exon_gff_df.iterrows():
        block_sizes = []
        chrom_starts = []
        exon_ids = []

        for i in range(len(row.ID_ex)):
            if pd.isna(row.ID_ex[i]):
                continue

            chrom_starts.append(str(row.start_ex[i]))
            block_sizes.append(str(row.end_ex[i] - row.start_ex[i]))
            exon_ids.append(
                (
                    row.exon_number[i],
                    orf_exon_map[
                        (row.start_ex[i],
                         row.end_ex[i],
                         row.strand_ex[i],
                         assembly_id)
                    ]
                )
            )

        block_sizes = ';'.join(block_sizes)
        chrom_starts = ';'.join(chrom_starts)
        transcript_idx_str = (
            f'{start}_{end}_{strand}_{assembly_id}_'
            f'{block_sizes}_{chrom_starts}'
        )
        transcript_idx = hashlib.sha3_512(
            transcript_idx_str.encode('utf-8')
        ).hexdigest()

        if transcript_idx in orf_transcript_idx_map.keys():
            transcript_orf_id = orf_transcript_idx_map[transcript_idx]

            update_entries.append({
                "id": transcript_orf_id,
                "chess_id": transcript_id,
                "attrs": {
                    "gene_type": row.gene_type_tx,
                    "db_xref": row.db_xref_tx,
                    "num_samples": row.num_samples_tx,
                    "max_tpm": row.max_tpm_tx,
                }
            })

    session.bulk_update_mappings(Transcript, update_entries)
    session.commit()


def update_transcript_utr_mapping(gencode_dir, version, session):
    """Update transcript UTR mappings.

    Args:
        gencode_dir (Path): Directory containing GENCODE files
        version (str): GENCODE version
        session: SQLAlchemy session object
    """
    gencode_expanded_gff = gencode_dir.joinpath(
        version,
        f'gencode.{version}.chr_patch_hapl_scaff.annotation.expanded.gff3'
    )

    gff_df = pd.read_csv(gencode_expanded_gff, sep='\t', low_memory=False)
    gff_df.fillna('', inplace=True)
    assembly_ids = {}

    for assembly in session.query(Assembly).all():
        if assembly.ucsc_style_name.startswith('chr') and \
           len(assembly.ucsc_style_name) < 6:
            assembly_ids[assembly.ucsc_style_name] = assembly.id
        else:
            assembly_ids[assembly.genbank_accession] = assembly.id

    gff_df['assembly_id'] = gff_df.apply(
        lambda x: assembly_ids[x.seq_id],
        axis=1
    )

    utr_five_gff_df = gff_df[
        gff_df['type'].isin(['five_prime_UTR'])
    ].copy()
    utr_three_gff_df = gff_df[
        gff_df['type'].isin(['three_prime_UTR'])
    ].copy()

    grouped_utr_five_df = utr_five_gff_df.groupby(
        ['start', 'end', 'strand', 'assembly_id']
    ).aggregate(list)
    grouped_utr_three_df = utr_three_gff_df.groupby(
        ['start', 'end', 'strand', 'assembly_id']
    ).aggregate(list)

    orf_ensembl_tx_map = {
        t.xref: t.transcript_id
        for t in session.query(TranscriptXref).filter(
            TranscriptXref.xref_dataset_id == 1
        ).all()
    }
    orf_utr5_idx_map = {
        (u.start, u.end, u.strand, u.assembly_id): u.id
        for u in session.query(UtrFive).all()
    }
    orf_utr3_idx_map = {
        (u.start, u.end, u.strand, u.assembly_id): u.id
        for u in session.query(UtrThree).all()
    }

    utr_five_txs = []
    utr_tx_entries = set()

    for idx, row in grouped_utr_five_df.iterrows():
        utr_orf_id = orf_utr5_idx_map[idx]
        for i, transcript_id in enumerate(row.transcript_id):
            transcript_orf_id = orf_ensembl_tx_map[transcript_id]
            
            if not (utr_orf_id, transcript_orf_id) in utr_tx_entries:
                utr_tx_entries.add((utr_orf_id, transcript_orf_id))
                utr_five_txs.append(
                    TranscriptUtrFive(
                        transcript_orf_id,
                        utr_orf_id,
                        attrs={
                            "tag": row.tag[i],
                            "transcript_name": row.transcript_name[i],
                            "exon_number": row.exon_number[i],
                            "transcript_support_level":
                                row.transcript_support_level[i]
                        }
                    )
                )

    session.add_all(utr_five_txs)
    session.commit()

    utr_three_txs = []
    utr_tx_entries = set()

    for idx, row in grouped_utr_three_df.iterrows():
        utr_orf_id = orf_utr3_idx_map[idx]
        for i, transcript_id in enumerate(row.transcript_id):
            transcript_orf_id = orf_ensembl_tx_map[transcript_id]
            
            if not (utr_orf_id, transcript_orf_id) in utr_tx_entries:
                utr_tx_entries.add((utr_orf_id, transcript_orf_id))
                utr_three_txs.append(
                    TranscriptUtrThree(
                        transcript_orf_id,
                        utr_orf_id,
                        attrs={
                            "tag": row.tag[i],
                            "transcript_name": row.transcript_name[i],
                            "exon_number": row.exon_number[i],
                            "transcript_support_level":
                                row.transcript_support_level[i]
                        }
                    )
                )
                
    session.add_all(utr_three_txs)
    session.commit()