"""Module to load core annotation data.

This module provides functions for loading various types of annotation data
into the ORF database, including genome assemblies, genes, transcripts,
exons, UTRs, and CDS regions from different data sources.
"""

from Bio.Seq import Seq
from orfdb.base import *  # Changed from veliadb.base
from orfdb import util as orf_utils  # Changed from veliadb import util

import ast
import gzip
import hashlib
import logging
import pandas as pd
from os.path import basename
from warnings import warn


def load_genome_assembly(session, assembly_df, genome_accession):
    """Load genome assembly information into database.

    Args:
        session: SQLAlchemy session object
        assembly_df (pd.DataFrame): Genome assembly report dataframe
        genome_accession (str): Genome assembly accession ID
    """
    assemblies = []
    for i, row in assembly_df.iterrows():
        assemblies.append(Assembly(
            row['GenBank-Accn'],
            row['RefSeq-Accn'],
            row['UCSC-style-name'],
            row['Sequence-Role'],
            row['Assembly-Unit'],
            row['Assigned-Molecule'],
            row['Assigned-Molecule-Location/Type'],
            row['Sequence-Length'],
            genome_accession=genome_accession
        ))
    
    # Add unmapped assembly entry
    assemblies.append(Assembly(
        genbank_accession="unmapped",
        refseq_accession="unmapped",
        ucsc_style_name="unmapped",
        sequence_role="unmapped",
        assembly_unit="unmapped",
        assigned_molecule="unmapped",
        assigned_molecule_location="unmapped",
        sequence_length=0,
        genome_accession="unmapped"
    ))

    session.add_all(assemblies)
    session.commit()


def load_gencode_genes(session, gene_gff_df, assembly_ids):
    """Load genes from GENCODE GFF file.

    Args:
        session: SQLAlchemy session object
        gene_gff_df (pd.DataFrame): GENCODE GFF dataframe filtered for genes
        assembly_ids (dict): Mapping of sequence IDs to assembly IDs
    """
    gene_gff_df['assembly_id'] = gene_gff_df.apply(
        lambda x: assembly_ids[x.seq_id],
        axis=1
    )

    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    gene_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_gene_gff_df = gene_gff_df.groupby(unique_cols).aggregate(list)

    genes = []
    synonym_dict = {}

    for idx, row in grouped_gene_gff_df.iterrows():
        genes.append(Gene(
            idx[0], idx[1], idx[2], idx[3],
            row.hgnc_id[0], row.gene_name[0], row.ID[0],
            '', '', '', row.gene_type[0], '', {}
        ))
        
        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(ID)

        syn_dict['HGNC_ID'] = set(row.hgnc_id)
        syn_dict['HGNC_SYMBOL'] = set(row.gene_name)
        synonym_dict[idx] = syn_dict

    # Add unmapped gene entry
    unmapped_assembly = session.query(Assembly).filter(
        Assembly.genbank_accession == 'unmapped'
    ).one()
    genes.append(Gene(
        start=-1, end=-1, strand='',
        assembly_id=unmapped_assembly.id,
        hgnc_id="unmapped",
        hgnc_name="unmapped",
        ensembl_id="unmapped",
        refseq_id="unmapped",
        velia_id="unmapped",
        chess_id="unmapped",
        gene_type="unmapped"
    ))

    session.add_all(genes)
    session.commit()

    logging.info(f'Added {len(genes)} GENCODE genes')

    # Add gene cross-references
    gene_xrefs = []
    ensembl_vdb_gene_map = {
        (g.start, g.end, g.strand, g.assembly_id): g.id
        for g in session.query(Gene).all()
    }
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '':
                    continue
                gene_xrefs.append(SequenceRegionXref(
                    ensembl_vdb_gene_map[idx],
                    synonym,
                    'synonym',
                    dataset_ids['ENSEMBL'],
                    dataset_ids[dataset_name]
                ))

    session.add_all(gene_xrefs)
    session.commit()

    logging.info(f'Added {len(gene_xrefs)} GENCODE gene synonyms')


def load_gencode_exons(session, exon_gff_df, assembly_ids):
    """Load exons from GENCODE GFF file.

    Args:
        session: SQLAlchemy session object
        exon_gff_df (pd.DataFrame): GENCODE GFF dataframe filtered for exons
        assembly_ids (dict): Mapping of sequence IDs to assembly IDs
    """
    exon_gff_df['assembly_id'] = exon_gff_df.apply(
        lambda x: assembly_ids[x.seq_id],
        axis=1
    )

    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    exon_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_exon_gff_df = exon_gff_df.groupby(unique_cols).aggregate(list)
        
    exons = []
    synonym_dict = {}

    for idx, row in grouped_exon_gff_df.iterrows():
        exons.append(Exon(
            idx[0], idx[1], idx[2], idx[3],
            row.exon_id[0], '', '', ''
        ))

        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.exon_id[i])
            syn_dict[row.source[i]].add(ID)

        synonym_dict[idx] = syn_dict

    session.add_all(exons)
    session.commit()

    logging.info(f'Added {len(exons)} GENCODE exons')
    
    exon_xrefs = []
    ensembl_orf_exon_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id 
        for e in session.query(Exon).all()
    }
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}
    
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '':
                    continue
                exon_xrefs.append(SequenceRegionXref(
                    ensembl_orf_exon_map[idx],
                    synonym,
                    'synonym',
                    dataset_ids['ENSEMBL'],
                    dataset_ids[dataset_name]
                ))

    session.add_all(exon_xrefs)
    session.commit()

    logging.info(f'Added {len(exon_xrefs)} GENCODE exon synonyms')


def load_gencode_transcripts(session, tx_gff_df, exon_gff_df, gencode_dir,
                           version, assembly_ids):
    """Load transcripts from GENCODE GFF file.

    Args:
        session: SQLAlchemy session object
        tx_gff_df (pd.DataFrame): GENCODE GFF dataframe filtered for transcripts
        exon_gff_df (pd.DataFrame): GENCODE GFF dataframe filtered for exons
        gencode_dir (Path): Path to GENCODE directory
        version (str): GENCODE version
        assembly_ids (dict): Mapping of sequence IDs to assembly IDs

    Returns:
        dict: Mapping of transcript IDs to exon information
    """
    tx_gff_df['assembly_id'] = tx_gff_df.apply(
        lambda x: assembly_ids[x.seq_id],
        axis=1
    )
    exon_gff_df['assembly_id'] = exon_gff_df.apply(
        lambda x: assembly_ids[x.seq_id],
        axis=1
    )

    merged_df = tx_gff_df.merge(
        exon_gff_df,
        left_on=('ID', 'seq_id'),
        right_on=('transcript_id', 'seq_id'),
        suffixes=('_tx', '_ex'),
        how='left'
    )

    transcript_exon_gff_df = merged_df.groupby(
        ['ID_tx', 'start_tx', 'end_tx', 'strand_tx', 'assembly_id_tx']
    ).aggregate(list)

    ensembl_orf_exon_map = {
        (e.start, e.end, e.strand, e.assembly_id): e.id
        for e in session.query(Exon).all()
    }
    ensembl_orf_gene_map = {
        g.ensembl_id: g.id
        for g in session.query(Gene).all()
    }

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

    unmapped_gene = session.query(Gene).filter(
        Gene.hgnc_id == 'unmapped'
    ).one()

    for (transcript_id, start, end, strand, assembly_id), row in \
            transcript_exon_gff_df.iterrows():
        block_sizes = []
        chrom_starts = []
        exon_ids = []

        try:
            gene_id = ensembl_orf_gene_map[row.gene_id_tx[0]]
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
                row.exon_number_ex[i],
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
                row.havana_transcript_tx
            )
            synonym_dict[transcript_idx]['CCDS'].update(row.ccdsid_tx)
        else:
            existing_entries.append(transcript_idx)
            synonym_dict[transcript_idx] = {}
            synonym_dict[transcript_idx]['ENSEMBL'] = set([transcript_id])
            synonym_dict[transcript_idx]['HAVANA'] = set(
                row.havana_transcript_tx
            )
            synonym_dict[transcript_idx]['CCDS'] = set(row.ccdsid_tx)

            transcripts.append(Transcript(
                start, end, strand, assembly_id,
                block_sizes, chrom_starts,
                transcript_idx, transcript_idx_str,
                gene_id, transcript_id, refseq_id, '', '',
                row.transcript_support_level_tx[0],
                row.transcript_type_tx[0],
                attrs={'tag': row.tag_tx[0]}
            ))

            transcript_exon_map[transcript_idx] = []

        for exon_num, exon_id in exon_ids:
            transcript_exon_map[transcript_idx].append((exon_num, exon_id))

    session.add_all(transcripts)
    session.commit()

    logging.info(f'Added {len(transcripts)} GENCODE transcripts')

    transcript_xrefs = []
    ensembl_orf_tx_map = {
        t.transcript_idx: t.id
        for t in session.query(Transcript).all()
    }
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '':
                    continue
                transcript_xrefs.append(TranscriptXref(
                    ensembl_orf_tx_map[idx],
                    synonym,
                    'synonym',
                    dataset_ids['ENSEMBL'],
                    dataset_ids[dataset_name]
                ))

    session.add_all(transcript_xrefs)
    session.commit()

    logging.info(f'Added {len(transcript_xrefs)} GENCODE transcript xrefs')

    return transcript_exon_map


def load_gencode_utr(session, utr_gff_df, assembly_ids, utr_class):
    """Load utrs from GENCODE gff

    Args:
        session (sqlalchemy.Session): An open database connection object
        utr_gff_df (pandas.DataFrame): Subset of the gencode gff with utr records
        assembly_ids (dict): Mapping of GENCODE 'seq_id' to internal database ids
    """

    utr_gff_df['assembly_id'] = utr_gff_df.apply(lambda x: assembly_ids[x.seq_id], axis=1)

    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    utr_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_utr_gff_df = utr_gff_df.groupby(unique_cols).aggregate(list)
        
    utrs = []
    synonym_dict = {}

    for idx, row in grouped_utr_gff_df.iterrows():

        utrs.append(utr_class(idx[0], idx[1], idx[2], idx[3],
                              row.ID[0], '', '', ''))
        
        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(ID)

        synonym_dict[idx] = syn_dict


    session.add_all(utrs)
    session.commit()

    logging.info(f'Added {len(utrs)} GENCODE UTRs')
    
    utr_xrefs = []

    ensembl_vdb_utr_map = {(u.start, u.end, u.strand, u.assembly_id): u.id for u in session.query(utr_class).all()}
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}
    
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                utr_xrefs.append(SequenceRegionXref(ensembl_vdb_utr_map[idx], synonym, 'synonym', dataset_ids['ENSEMBL'], dataset_ids[dataset_name]))

    session.add_all(utr_xrefs)
    session.commit()

    logging.info(f'Added {len(utr_xrefs)} GENCODE UTR synonyms')


def load_gencode_cds(session, cds_gff_df, assembly_ids):
    """Load CDS from GENCODE gff

    Args:
        session (sqlalchemy.Session): An open database connection object
        cds_gff_df (pandas.DataFrame): Subset of the gencode gff with CDS records
        assembly_ids (dict): Mapping of GENCODE 'seq_id' to internal database ids
    """

    cds_gff_df['assembly_id'] = cds_gff_df.apply(lambda x: assembly_ids[x.seq_id], axis=1)

    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    cds_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_cds_gff_df = cds_gff_df.groupby(unique_cols).aggregate(list)
    
    cds = []
    synonym_dict = {}
    
    for idx, row in grouped_cds_gff_df.iterrows():
        
        synonym_dict[idx] = {'ENSEMBL': set()}
        
        for i in range(len(row.ID)):
            ensembl_id = f'cds-{row.transcript_id[i]}-{row.protein_id[i]}-{row.exon_number[i]}'
            
            if i == 0:
                cds.append(Cds(idx[0], idx[1], idx[2], idx[3],
                            ensembl_id, '', '', '', 
                            row.ccdsid[i], row.protein_id[i], ''))

            synonym_dict[idx]['ENSEMBL'].add(ensembl_id)
                
    session.add_all(cds)
    session.commit()

    logging.info(f'Added {len(cds)} GENCODE CDS')

    cds_xrefs = []

    ensembl_vdb_cds_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}
    
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                cds_xrefs.append(SequenceRegionXref(ensembl_vdb_cds_map[idx], synonym, 'synonym', dataset_ids['ENSEMBL'], dataset_ids[dataset_name]))

    session.add_all(cds_xrefs)
    session.commit()

    logging.info(f'Added {len(cds_xrefs)} GENCODE CDS synonyms')


def load_gencode_orfs(session, orf_gff_df, assembly_ids):
    """Load ORFs from GENCODE gff

    Args:
        session (sqlalchemy.Session): An open database connection object
        orf_gff_df (pandas.DataFrame): Subset of the gencode gff with CDS, start_codon, or stop_codon records
        assembly_ids (dict): Mapping of GENCODE 'seq_id' to internal database ids
    """

    orf_gff_df['assembly_id'] = orf_gff_df.apply(lambda x: assembly_ids[x.seq_id], axis=1)
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    orf_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_orf_gff_df = orf_gff_df.groupby('protein_id').aggregate(list)

    ensembl_vdb_cds_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}

    orfs = []
    cds_orf_map = {}
    synonym_dict = {}

    for protein_id, row in grouped_orf_gff_df.iterrows():
        block_sizes = []
        chrom_starts = []
        exon_frames = []
        cds_ids = []
        start = row.start[0]
        end = row.end[0]
        
        for i in range(len(row.ID)):
            chrom_starts.append(str(row.start[i]))
            block_sizes.append(str(row.end[i] - row.start[i]))
            exon_frames.append(str(row.phase[i]))
            
            if row.start[i] < start:
                start = row.start[i]
            if row.end[i] > end:
                end = row.end[i]
                
            cds_ids.append((row.exon_number[i], row.phase[i], ensembl_vdb_cds_map[(row.start[i], row.end[i], row.strand[i], row.assembly_id[i])]))
        
        block_sizes = ';'.join(block_sizes)
        chrom_starts = ';'.join(chrom_starts)
        exon_frames = ';'.join(exon_frames)
        orf_val = f'{start}_{end}_{row.strand[0]}_{row.assembly_id[0]}_{block_sizes}_{chrom_starts}_{exon_frames}'
        orf_idx = hashlib.sha3_512(orf_val.encode('utf-8')).hexdigest()

        if orf_idx in synonym_dict.keys():
            synonym_dict[orf_idx]['ENSEMBL'].add(protein_id)
            continue
        else:
            synonym_dict[orf_idx] = {'ENSEMBL': set([protein_id])}

        orfs.append(Orf(start, end, row.strand[0], row.assembly_id[0], 
                        block_sizes, chrom_starts, exon_frames, orf_idx, 
                        orf_val, '', '', '', protein_id, '', ''))

        cds_orf_map.update({cds_id: (orf_idx, exon_num, phase) \
                                        for (exon_num, phase, cds_id) in cds_ids})

    session.add_all(orfs)
    session.commit()

    logging.info(f'Added {len(orfs)} GENCODE ORFs')

    orf_xrefs = []

    ensembl_vdb_orf_map = {o.orf_idx: o.id for o in session.query(Orf).all()}

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}
    
    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '': continue
                orf_xrefs.append(OrfXref(ensembl_vdb_orf_map[idx], synonym, 'synonym', dataset_ids['ENSEMBL'], dataset_ids[dataset_name]))

    session.add_all(orf_xrefs)
    session.commit()

    logging.info(f'Added {len(orf_xrefs)} GENCODE ORF synonyms')
    
    return cds_orf_map


def load_gencode_orf_cds(session, cds_orf_map):
    """Load mappings of ORF <-> CDS from GENCODE gff

    Args:
        session (sqlalchemy.Session): An open database connection object
        orf_cds_map (dict): Mapping of orf unique_id (e.g. ensembl_protein_id) to CDSs
    """

    cds_orfs = []

    ensembl_orf_ids = {o.orf_idx: o.id for o in session.query(Orf).all()}

    for i, (cds_id, orf_idx) in enumerate(cds_orf_map.items()):
        orf_id = ensembl_orf_ids[orf_idx[0]]
        
        cds_orfs.append(CdsOrf(cds_id, orf_id, orf_idx[-2], orf_idx[-1]))
        
    session.add_all(cds_orfs)
    session.commit()


def load_gencode_exon_cds(session, cds_gff_df, assembly_ids):
    """
    """
    ensembl_vdb_cds_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}
    ensembl_vdb_exon_map = {(exon.assembly_id, syn.xref): exon.id for syn, exon in session.query(SequenceRegionXref, Exon)\
                            .join(Exon, Exon.id == SequenceRegionXref.sequence_region_id).all()}

    cds_gff_df['assembly_id'] = cds_gff_df.apply(lambda x: assembly_ids[x.seq_id], axis=1)

    exon_cds_ids = []
    exon_cds = []

    for i, row in cds_gff_df.iterrows():
        
        exon_cds_ids.append((ensembl_vdb_exon_map[(row.assembly_id, row.exon_id)], 
                             ensembl_vdb_cds_map[(row.start, row.end, row.strand, row.assembly_id)]))

    x = set(exon_cds_ids).difference([(ec.exon_id, ec.cds_id) for ec in session.query(ExonCds).all()])    

    for (exon_id, cds_id) in x:
        exon_cds.append(ExonCds(exon_id, cds_id))

    session.add_all(exon_cds)
    session.commit()


def load_gencode_riboseq_orfs(session, orf_bed_df, assembly_ids, genome_seq):
    """Load GENCODE auxiliary Ribo-seq ORFs (Mudge et. al, Nature Biotech)

    Args:
        session (sqlalchemy.Session): An open database connection object
        orf_bed_df (pd.DataFrame): Ribo-seq ORF bed file from GENCODE
        assembly_ids (dict): Mapping of GENCODE 'seq_id' to internal database ids

    """

    ensembl_vdb_transcript_map = {t.ensembl_id.split('.')[0]: t.id for t in session.query(Transcript).all()}
    ensembl_vdb_cds_map = {(c.start, c.end, c.strand, c.assembly_id): c.ensembl_id for c in session.query(Cds).all()}
    ensembl_vdb_exon_map = {(e.start, e.end, e.strand, e.assembly_id): e.ensembl_id for e in session.query(Exon).all()}

    orf_entries = []
    cds_entries = []
    exon_entries = []

    missed_transcripts = []

    transcript_exon_map = {}
    transcript_orf_map = {}
    exon_cds_map = {}
    cds_orf_map = {}

    for i, row in orf_bed_df.iterrows():

        transcript_ids = []

        assembly_id = assembly_ids[row.chrom]
        orf_val = f'{row.chromStart}_{row.chromEnd}_{row.strand}_{assembly_id}_{row.blockSizes}_{row.chromStarts}_{row.exonFrames}'
        orf_idx = hashlib.sha3_512(orf_val.encode('utf-8')).hexdigest()

        for tx_id in row.all_transcript_ids.split(','):
            try:
                transcript_ids.append(ensembl_vdb_transcript_map[tx_id])
            except:
                missed_transcripts.append(tx_id)

        nt_seqs = []
        chrom = genome_seq[row.chrom]
        exon_frames = row.exonFrames.rstrip(',').split(',')

        for i, (block, cstart) in enumerate(zip(row.blockSizes.rstrip(',').split(','), row.chromStarts.rstrip(',').split(',')), start=1):

            cds_start = int(row.chromStart) + int(cstart)
            cds_end = cds_start + int(block)
            
            cds_idx = (cds_start, cds_end, row.strand, assembly_ids[row.chrom])
            
            if cds_idx in ensembl_vdb_cds_map.keys():
                cds_id = ensembl_vdb_cds_map[cds_idx]
            else:
                cds_id = f"cds-{row['name']}-{i}"
                cds_entries.append(Cds(cds_start, cds_end, row.strand, assembly_ids[row.chrom], cds_id, '', '', '', '', '', ''))
                ensembl_vdb_cds_map[cds_idx] = cds_id

            if cds_idx in ensembl_vdb_exon_map.keys():
                exon_id = ensembl_vdb_exon_map[cds_idx]
            else:
                exon_id = f"exon-{row['name']}-{i}"
                exon_entries.append(Exon(cds_start, cds_end, row.strand, assembly_ids[row.chrom], exon_id, '', '', ''))
                ensembl_vdb_exon_map[cds_idx] = exon_id
                
            for transcript_id in transcript_ids:
                if transcript_id in transcript_exon_map.keys():
                    transcript_exon_map[transcript_id].append((exon_id, i))
                else:
                    transcript_exon_map[transcript_id] = [(exon_id, i)]

            exon_cds_map[exon_id] = cds_id
            cds_orf_map[cds_id] = (orf_idx, i, exon_frames[i-1])

            if row.strand == '+':
                nt_seqs.append(str(chrom[cds_start:cds_end].seq))

            elif row.strand == '-':
                nt_seqs.insert(0, str(chrom[cds_start:cds_end].reverse_complement().seq))

        assert row.sequence == str(Seq(''.join(nt_seqs)).translate()[:-1])

        for transcript_id in transcript_ids:
            transcript_orf_map[transcript_id] = orf_idx

        orf_entries.append(Orf(row.chromStart, row.chromEnd, row.strand, assembly_ids[row.chrom],
                               row.blockSizes, row.chromStarts, row.exonFrames, orf_idx, orf_val, row['name'], '', str(Seq(''.join(nt_seqs))), '', '', ''))


    session.add_all(cds_entries)
    session.commit()

    session.add_all(exon_entries)
    session.commit()

    session.add_all(orf_entries)
    session.commit()

    return (transcript_exon_map, transcript_orf_map, exon_cds_map, cds_orf_map)


def load_gencode_riboseq_orf_cds(session, cds_orf_map):
    """
    """
    ensembl_vdb_cds_map = {c.ensembl_id: c.id for c in session.query(Cds).all()}
    ensembl_vdb_orf_map = {o.orf_idx: o.id for o in session.query(Orf).all()}

    cds_orfs = []

    for (cds_id, orf_idx) in cds_orf_map.items():
        cds_orfs.append(CdsOrf(ensembl_vdb_cds_map[cds_id], ensembl_vdb_orf_map[orf_idx[0]], orf_idx[-2], orf_idx[-1]))

    session.add_all(cds_orfs)
    session.commit()


def load_gencode_riboseq_exon_cds(session, exon_cds_map):
    """
    """
    ensembl_vdb_cds_map = {c.ensembl_id: c.id for c in session.query(Cds).all()}
    ensembl_vdb_exon_map = {e.ensembl_id: e.id for e in session.query(Exon).all()}

    exon_cds = []

    for exon, cds in exon_cds_map.items():
        exon_cds.append(ExonCds(ensembl_vdb_exon_map[exon], ensembl_vdb_cds_map[cds]))

    session.add_all(exon_cds)
    session.commit()


def load_gencode_riboseq_transcript_exons(session, transcript_exon_rb_map):
    """
    """
    transcript_exon_entries = [(te.transcript_id, te.exon_id, te.exon_number) for te in session.query(TranscriptExon).all()]
    ensembl_vdb_exon_map = {e.ensembl_id: e.id for e in session.query(Exon).all()}

    transcript_exons = []
    all_entries = []

    for transcript_id, vals in transcript_exon_rb_map.items():
        for (exon_id, exon_number) in vals:
            all_entries.append((transcript_id, ensembl_vdb_exon_map[exon_id], exon_number))
            
    unique_entries = set(all_entries).difference(set(transcript_exon_entries))

    for entry in unique_entries:
        transcript_exons.append(TranscriptExon(*entry))
        
    session.add_all(transcript_exons)
    session.commit()


def load_gencode_riboseq_transcript_orfs(session, transcript_orf_map):
    """
    """
    ensembl_vdb_orf_map = {o.orf_idx: o.id for o in session.query(Orf).all()}

    transcript_orfs = []

    for transcript_id, orf_idx in transcript_orf_map.items():
        transcript_orfs.append(TranscriptOrf(transcript_id, ensembl_vdb_orf_map[orf_idx], 'GENCODE_riboseq'))

    session.add_all(transcript_orfs)
    session.commit()


def load_gencode_transcript_exons(session, transcript_exon_map):
    """Load transcript-exon relationships from GENCODE data.

    Args:
        session: SQLAlchemy session object
        transcript_exon_map (dict): Mapping of transcript IDs to exon info
    """
    transcript_exons = []
    ensembl_orf_tx_map = {
        t.transcript_idx: t.id
        for t in session.query(Transcript).all()
    }

    for transcript_idx, exon_entries in transcript_exon_map.items():
        for exon_num, exon_id in exon_entries:
            transcript_exons.append(TranscriptExon(
                ensembl_orf_tx_map[transcript_idx],
                exon_id,
                exon_num
            ))

    session.add_all(transcript_exons)
    session.commit()

    logging.info(f'Added {len(transcript_exons)} GENCODE transcript exons')


def load_gencode_transcript_orfs(session, cds_gff_df):
    """Load transcript-ORF relationships from GENCODE data.

    Args:
        session: SQLAlchemy session object
        cds_gff_df (pd.DataFrame): GENCODE GFF dataframe filtered for CDS
    """
    transcript_orfs = []
    ensembl_orf_tx_map = {
        t.ensembl_id: t.id
        for t in session.query(Transcript).all()
    }
    ensembl_orf_orf_map = {
        o.ensembl_id: o.id
        for o in session.query(Orf).all()
    }

    for i, row in cds_gff_df.iterrows():
        if row.protein_id not in ensembl_orf_orf_map.keys():
            continue
        if row.transcript_id not in ensembl_orf_tx_map.keys():
            continue

        transcript_orfs.append(TranscriptOrf(
            ensembl_orf_tx_map[row.transcript_id],
            ensembl_orf_orf_map[row.protein_id]
        ))

    session.add_all(transcript_orfs)
    session.commit()

    logging.info(f'Added {len(transcript_orfs)} GENCODE transcript ORFs')


def load_gencode_orf_cds(session, cds_orf_map):
    """Load ORF-CDS relationships from GENCODE data.

    Args:
        session: SQLAlchemy session object
        cds_orf_map (dict): Mapping of ORF IDs to CDS information
    """
    orf_cds = []
    ensembl_orf_orf_map = {
        o.orf_idx: o.id
        for o in session.query(Orf).all()
    }

    for orf_idx, cds_entries in cds_orf_map.items():
        for exon_num, phase, cds_id in cds_entries:
            orf_cds.append(OrfCds(
                ensembl_orf_orf_map[orf_idx],
                cds_id,
                exon_num,
                phase
            ))

    session.add_all(orf_cds)
    session.commit()

    logging.info(f'Added {len(orf_cds)} GENCODE ORF CDS')


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


def load_refseq_genes(session, gene_gff_df):
    """
    """

    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    gene_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_gene_gff_df = gene_gff_df.groupby(unique_cols).aggregate(list)

    vdb_gene_idx_map = {(g.start, g.end, g.strand, g.assembly_id): g.id for g in session.query(Gene).all()}

    update_entries = []
    refseq_only_genes = []

    synonym_dict = {}

    for gene_idx, row in grouped_gene_gff_df.iterrows():

        if gene_idx in vdb_gene_idx_map.keys():
            gene_id = vdb_gene_idx_map[gene_idx]

            update_entries.append({"id": gene_id, "refseq_id": row.entrez_gene_id[0],
                                   "attrs": {"description": row.description[0],
                                            "gene_biotype": row.gene_biotype[0]}})

        elif gene_idx not in synonym_dict.keys():
            refseq_only_genes.append(Gene(gene_idx[0], gene_idx[1], gene_idx[2], gene_idx[3], 
                                        row.HGNC_ID[0], row.Name[0], '', row.entrez_gene_id[0], '', '', row.gene_biotype[0], 
                                        attrs={'description': row.description[0]}))
                    
        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.Name[i])
            syn_dict[row.source[i]].add(ID)
            syn_dict[row.source[i]].add(row.entrez_gene_id[i])

        syn_dict['HGNC_ID'] = set(row.HGNC_ID)
        synonym_dict[gene_idx] = syn_dict

    session.bulk_update_mappings(Gene, update_entries,)
    session.commit()

    session.add_all(refseq_only_genes)
    session.commit()

    logging.info(f'Updated {len(update_entries)} genes with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(refseq_only_genes)} RefSeq genes without an exact coordinate match to GENCODE')

    gene_xrefs = []

    vdb_gene_idx_map = {(g.start, g.end, g.strand, g.assembly_id): g.id for g in session.query(Gene).all()}    
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '': continue
                gene_xrefs.append(SequenceRegionXref(vdb_gene_idx_map[idx], synonym, 'synonym', dataset_ids['RefSeq'], dataset_ids[dataset_name]))

    session.add_all(gene_xrefs)
    session.commit()

    logging.info(f'Added {len(gene_xrefs)} RefSeq gene synonyms')


def load_refseq_transcripts(session, tx_gff_df, exon_gff_df):
    """
    """
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
                'inference', 'Note', 'experiment', 'gbkey',] 
        
    merged_df = tx_gff_df[map_cols].merge(exon_gff_df[map_cols + ['exon_number']], 
                                          left_on=('ID', 'seq_id'), 
                                          right_on=('Parent', 'seq_id'),
                                          suffixes=('_tx', '_ex'), how='left')    

    merged_df.sort_values(by=['ID_tx', 'start_ex'], inplace=True)

    transcript_exon_gff_df = merged_df.groupby(['ID_tx', 'start_tx', 'end_tx', 'strand_tx', 'assembly_id_tx']).aggregate(list)

    vdb_exon_map = {(e.start, e.end, e.strand, e.assembly_id): e.id for e in session.query(Exon).all()}
    refseq_vdb_gene_map = {g.refseq_id: g.id for g in session.query(Gene).all()}
    vdb_transcript_idx_map = {t.transcript_idx: t.id for t in session.query(Transcript).all()}

    transcripts = []
    transcript_exon_map = {}

    synonym_dict = {}

    update_entries = []
    existing_entries = []

    unmapped_gene = session.query(Gene).filter(Gene.hgnc_id == 'unmapped').one()

    for (transcript_id, start, end, strand, assembly_id), row in transcript_exon_gff_df.iterrows():
        block_sizes = []
        chrom_starts = []
        exon_ids = []

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

            update_entries.append({"id": transcript_vdb_id, "refseq_id": transcript_id,
                                   "attrs": {"inference": row.inference_tx[0],
                                             "experiment": row.experiment_tx[0],
                                             "product": row['product_tx'][0],
                                             "parent": row.Parent_tx[0],
                                             "gbkey": row.gbkey_tx[0]}})

        elif transcript_idx not in existing_entries:
            existing_entries.append(transcript_idx)
            
            transcripts.append(Transcript(start, end, strand, assembly_id,
                                        block_sizes, chrom_starts, transcript_idx, transcript_idx_str,
                                        gene_id, '', transcript_id, '', '', 
                                        row.model_evidence_tx[0], row.type_tx[0], 
                                        attrs={"inference": row.inference_tx[0],
                                               "experiment": row.experiment_tx[0],
                                               "product": row['product_tx'][0],
                                               "parent": row.Parent_tx[0],
                                               "gbkey": row.gbkey_tx[0]}))
        
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

    session.bulk_update_mappings(Transcript, update_entries,)
    session.commit()

    session.add_all(transcripts)
    session.commit()

    logging.info(f'Updated {len(update_entries)} transcripts with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(transcripts)} RefSeq transcripts without an exact coordinate match to GENCODE')

    transcript_xrefs = []

    ensembl_vdb_tx_map = {t.transcript_idx: t.id for t in session.query(Transcript).all()}

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '': continue
                transcript_xrefs.append(TranscriptXref(ensembl_vdb_tx_map[idx], synonym, 'synonym', dataset_ids['RefSeq'], dataset_ids[dataset_name]))

    session.add_all(transcript_xrefs)
    session.commit()

    return transcript_exon_map


def load_refseq_exons(session, exon_gff_df):
    """
    """
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    exon_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_exon_gff_df = exon_gff_df.groupby(unique_cols).aggregate(list)

    vdb_exon_idx_map = {(e.start, e.end, e.strand, e.assembly_id): e.id for e in session.query(Exon).all()}

    exons = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    for exon_idx, row in grouped_exon_gff_df.iterrows():
        
        if exon_idx in vdb_exon_idx_map.keys():
            exon_id = vdb_exon_idx_map[exon_idx]
            
            update_entries.append({"id": exon_id, "refseq_id": row.ID[0],
                                "attrs": {"product": row['product'][0],
                                            "parent": row.Parent[0]}})
            
        elif exon_idx not in existing_entries:
            
            existing_entries.append(exon_idx)
        
            exons.append(Exon(exon_idx[0], exon_idx[1], exon_idx[2], exon_idx[3],
                            '', row.ID[0], '', '', 
                            attrs={"inference": row.inference[0],
                                    "model_evidence": row.model_evidence[0],
                                    "notes": row.Note[0]}))
        
        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.Name[i])
            syn_dict[row.source[i]].add(ID)

        synonym_dict[exon_idx] = syn_dict

    session.bulk_update_mappings(Exon, update_entries,)
    session.commit()

    session.add_all(exons)
    session.commit()

    logging.info(f'Updated {len(update_entries)} exons with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(exons)} RefSeq exon without an exact coordinate match to GENCODE')

    vdb_exon_idx_map = {(e.start, e.end, e.strand, e.assembly_id): e.id for e in session.query(Exon).all()}

    exon_xrefs = []

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '': continue
                exon_xrefs.append(SequenceRegionXref(vdb_exon_idx_map[idx], synonym, 'synonym', dataset_ids['RefSeq'], dataset_ids[dataset_name]))

    session.add_all(exon_xrefs)
    session.commit()


def load_refseq_cds(session, cds_gff_df):
    """
    """
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    cds_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_cds_gff_df = cds_gff_df.groupby(unique_cols).aggregate(list)

    vdb_cds_idx_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}

    cds = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    for cds_idx, row in grouped_cds_gff_df.iterrows():
        
        if cds_idx in vdb_cds_idx_map.keys():
            cds_id = vdb_cds_idx_map[cds_idx]
            
            update_entries.append({"id": cds_id, "refseq_id": row.ID[0], "refseq_protein_id": row.protein_id[0]})
            
        elif cds_idx not in existing_entries:
            
            existing_entries.append(cds_idx)
        
            cds.append(Cds(cds_idx[0], cds_idx[1], cds_idx[2], cds_idx[3],
                        '', row.Name[0], '', '', '', '', row.protein_id[0]))

        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.Name[i])
            syn_dict[row.source[i]].add(ID)

        synonym_dict[cds_idx] = syn_dict

    session.bulk_update_mappings(Cds, update_entries,)
    session.commit()
    
    session.add_all(cds)
    session.commit()

    logging.info(f'Updated {len(update_entries)} CDS with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(cds)} RefSeq CDS without an exact coordinate match to GENCODE')

    vdb_cds_idx_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}

    cds_xrefs = []

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '': continue
                cds_xrefs.append(SequenceRegionXref(vdb_cds_idx_map[idx], synonym, 'synonym', dataset_ids['RefSeq'], dataset_ids[dataset_name]))

    session.add_all(cds_xrefs)
    session.commit()


def load_refseq_orfs(session, cds_gff_df):
    """
    """
    orf_gff_df = cds_gff_df[cds_gff_df['protein_id'] != ''].copy()
    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    orf_gff_df.sort_values(by=unique_cols, inplace=True)    
    orf_gff_df = orf_gff_df.groupby('protein_id').aggregate(list)

    refseq_vdb_cds_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}
    vdb_orf_idx_map = {o.orf_idx: o.id for o in session.query(Orf).all()}

    orfs = []
    cds_orf_map = {}
    existing_entries = []
    update_entries = []

    synonym_dict = {}

    for protein_id, row in orf_gff_df.iterrows():
        block_sizes = []
        chrom_starts = []
        exon_frames = []
        cds_ids = []
        start = row.start[0]
        end = row.end[0]

        for i in range(len(row.ID)):
            chrom_starts.append(str(row.start[i]))
            block_sizes.append(str(row.end[i] - row.start[i]))
            exon_frames.append(str(row.phase[i]))

            if row.start[i] < start:
                start = row.start[i]
            if row.end[i] > end:
                end = row.end[i]

            # TODO this is NOT necessarily the correct CDS number in terms of mapping to exon number
            # Since RefSeq does NOT provide this mapping, at this point in the code a position based
            # lookup would need to be referenced to get the correct exon number
            cds_ids.append((i, row.phase[i], refseq_vdb_cds_map[(row.start[i], row.end[i], row.strand[i], row.assembly_id[i])]))

        block_sizes = ';'.join(block_sizes)
        chrom_starts = ';'.join(chrom_starts)
        exon_frames = ';'.join(exon_frames)
        orf_val = f'{start}_{end}_{row.strand[0]}_{row.assembly_id[0]}_{block_sizes}_{chrom_starts}_{exon_frames}'
        orf_idx = hashlib.sha3_512(orf_val.encode('utf-8')).hexdigest()
        
        if orf_idx in vdb_orf_idx_map.keys():
            orf_id = vdb_orf_idx_map[orf_idx]
            
            update_entries.append({"id": orf_id, "refseq_protein_id": row.name})
        
        elif orf_idx not in existing_entries:
            
            existing_entries.append(orf_idx)

            orfs.append(Orf(start, end, row.strand[0], row.assembly_id[0], 
                        block_sizes, chrom_starts, exon_frames, orf_idx, orf_val,
                        '', '', '', protein_id, '', ''))
        
        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.Name[i])
            syn_dict[row.source[i]].add(ID)
            syn_dict[row.source[i]].add(protein_id)

        synonym_dict[orf_idx] = syn_dict
        
        cds_orf_map.update({cds_id: (orf_idx, exon_num, phase) \
                                        for (exon_num, phase, cds_id) in cds_ids})


    session.bulk_update_mappings(Orf, update_entries,)
    session.commit()

    session.add_all(orfs)
    session.commit()

    logging.info(f'Updated {len(update_entries)} ORFs with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(orfs)} RefSeq ORFs without an exact coordinate match to GENCODE')

    orf_xrefs = []

    vdb_orf_idx_map = {o.orf_idx: o.id for o in session.query(Orf).all()}

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '': continue
                orf_xrefs.append(OrfXref(vdb_orf_idx_map[idx], synonym, 'synonym', dataset_ids['RefSeq'], dataset_ids[dataset_name]))

    session.add_all(orf_xrefs)
    session.commit()

    return cds_orf_map


def load_refseq_orf_cds(session, cds_orf_map):
    """
    """
    cds_orfs = []

    vdb_orf_idx_map = {o.orf_idx: o.id for o in session.query(Orf).all()}
    vdb_cds_orf_map = {(co.cds_id, co.orf_id) for co in session.query(CdsOrf).all()}

    for i, (cds_id, orf_idx) in enumerate(cds_orf_map.items()):
        orf_id = vdb_orf_idx_map[orf_idx[0]]
        
        # TODO at some point we might trust RefSeq CDS numbers and may
        # want to enable update or at least tracking of the conflict 
        if (cds_id, orf_id) not in vdb_cds_orf_map:
            cds_orfs.append(CdsOrf(cds_id, orf_id, orf_idx[-2], orf_idx[-1]))

    session.add_all(cds_orfs)
    session.commit()


def load_refseq_lncRNAs(session, lnc_gff_df):
    """
    This function is almost identical to load transcripts.  In the GENCODE loading
    there is no distinction for transcripts/lncRNAs.  However, in the RefSeq GFF
    there are a few different fields and cleary a different annotation process that
    make lncRNAs look more like genes.  I still load them here as transcripts but
    it should be possible to merge this code with transcript loading in the future.
    """

    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    lnc_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_lnc_gff_df = lnc_gff_df.groupby(unique_cols).aggregate(list)

    synonym_vdb_gene_map = {gs.xref: gs.sequence_region_id for gs in session.query(SequenceRegionXref).join(Gene, Gene.id == SequenceRegionXref.sequence_region_id).all()}
    vdb_tx_idx_map = {(t.start, t.end, t.strand, t.assembly_id): t.id for t in session.query(Transcript).all()}

    transcripts = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    unmapped_gene = session.query(Gene).filter(Gene.hgnc_id == 'unmapped').one()

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

    session.bulk_update_mappings(Transcript, update_entries,)
    session.commit()

    session.add_all(transcripts)
    session.commit()

    logging.info(f'Updated {len(update_entries)} lncRNAs with RefSeq info that had exact coordinate matches to GENCODE')
    logging.info(f'Added {len(transcripts)} RefSeq lncRNAs without an exact coordinate match to GENCODE')

    vdb_seq_syn_map = {(ss.sequence_region_id, ss.xref, ss.xref_dataset_id) for ss in session.query(SequenceRegionXref).all()}
    vdb_tx_idx_map = {(t.start, t.end, t.strand, t.assembly_id): t.id for t in session.query(Transcript).all()}

    transcript_xrefs = []

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if (vdb_tx_idx_map[idx], synonym, dataset_ids['RefSeq']) not in vdb_seq_syn_map and synonym != '':
                    transcript_xrefs.append(SequenceRegionXref(vdb_tx_idx_map[idx], synonym, 'synonym', dataset_ids['RefSeq'], dataset_ids[dataset_name]))

    session.add_all(transcript_xrefs)
    session.commit()


def load_refseq_transcript_exons(session, transcript_exon_map):
    """
    """
    refseq_dataset_id = session.query(Dataset.id).filter(Dataset.name == 'RefSeq').one()[0]

    refseq_vdb_transcript_map = {t.refseq_id: t.id for t in session.query(Transcript).all()}

    refseq_vdb_transcript_map.update({tx.xref: tx.transcript_id for tx in \
                                     session.query(TranscriptXref)\
                                            .filter(and_(TranscriptXref.transcript_dataset_id == refseq_dataset_id,
                                                         TranscriptXref.xref_dataset_id == refseq_dataset_id)).all()})
    
    transcript_exons = []

    existing_entries = set([(te.transcript_id, te.exon_id, te.exon_number) for te in session.query(TranscriptExon).all()])
    entries = set()

    for transcript_id, vals in transcript_exon_map.items():
        for (exon_number, exon_id) in vals:
            entries.add((refseq_vdb_transcript_map[transcript_id], exon_id, exon_number))

    entries.difference_update(existing_entries)
            
    for entry in entries:
        transcript_exons.append(TranscriptExon(*entry, attrs={'mapping_source': 'RefSeq'}))
        
    session.add_all(transcript_exons)
    session.commit()

    logging.info(f'Added {len(transcript_exons)} transcript <-> exon mappings from RefSeq')


def load_refseq_transcript_orfs(session, cds_gff_df):
    """
    """
    vdb_orf_synonym_map = {o.xref: o.orf_id for o in session.query(OrfXref).all()}

    vdb_transcript_synonym_map = {(t.assembly_id, ts.xref): t.id for ts,t in session.query(SequenceRegionXref, Transcript)\
                                    .join(Transcript, Transcript.id == SequenceRegionXref.sequence_region_id).all()}

    vdb_tx_orf_ids = {(to.transcript_id, to.orf_id) for to in session.query(TranscriptOrf).all()}

    transcript_orfs = []
    missed_transcripts = []

    for i, row in cds_gff_df.iterrows():

        if row.protein_id == '' or row.Parent == '':
            continue
        
        refseq_tx_id = row.Parent
        
        try:
            orf_id = vdb_orf_synonym_map[row.protein_id]
            transcript_id = vdb_transcript_synonym_map[(row.assembly_id, refseq_tx_id)]

        except:
            missed_transcripts.append((row.ID, row.Parent))
            continue
            
        if (transcript_id, orf_id) not in vdb_tx_orf_ids:
            
            transcript_orfs.append(TranscriptOrf(transcript_id, orf_id, row['product']))
            vdb_tx_orf_ids.add((transcript_id, orf_id))

    session.add_all(transcript_orfs)
    session.commit()
    
    logging.info(f'Added {len(transcript_orfs)} transcript <-> orf mappings from RefSeq')
    logging.info(f'Missed {len(missed_transcripts)} potential mappings without a mappable transcript from the cds entry')


def load_refseq_exon_cds():
    """
    TODO: Additional code needs to be written to try and map
    RefSeq exons to RefSeq CDS.  For some reason this mapping
    is not provided in the GFF and so a position based lookup
    would be needed to enable loading this for RefSeq
    """
    return


def load_chess_exons(session, exon_gff_df, assembly_ids):
    """Load exon data from CHESS database.

    This function processes and loads exon annotations from the CHESS database,
    creating new exon entries and their corresponding cross-references.

    Args:
        session: SQLAlchemy session object
        exon_gff_df (pd.DataFrame): CHESS GFF dataframe filtered for exons.
            Expected columns include:
            - start: Exon start position
            - end: Exon end position
            - strand: Strand orientation
            - seq_id: Chromosome/sequence identifier
            - ID: Exon identifier
            - source: Data source
        assembly_ids (dict): Mapping of sequence IDs to assembly IDs

    Returns:
        None
    """

    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    exon_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_exon_gff_df = exon_gff_df.groupby(unique_cols).aggregate(list)

    vdb_exon_idx_map = {(e.start, e.end, e.strand, e.assembly_id): e.id for e in session.query(Exon).all()}

    exons = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    for exon_idx, row in grouped_exon_gff_df.iterrows():
        chess_id = f'exon-{row.Parent[0]}'

        if exon_idx in vdb_exon_idx_map.keys():
            exon_id = vdb_exon_idx_map[exon_idx]
            update_entries.append({"id": exon_id, "chess_id": chess_id,
                                   "attrs": {"parent": row.Parent[0]}})

        elif exon_idx not in existing_entries:

            existing_entries.append(exon_idx)
            exons.append(Exon(exon_idx[0], exon_idx[1], exon_idx[2], exon_idx[3],
                            '', '', chess_id, '', attrs={"parent": row.Parent[0]}))

        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(chess_id)

        synonym_dict[exon_idx] = syn_dict

    session.bulk_update_mappings(Exon, update_entries,)
    session.commit()

    session.add_all(exons)
    session.commit()

    logging.info(f'Updated {len(update_entries)} exons with CHESS info that had exact coordinate matches to GENCODE/RefSeq')
    logging.info(f'Added {len(exons)} CHESS exon without an exact coordinate match to GENCODE/RefSeq')

    vdb_exon_idx_map = {(e.start, e.end, e.strand, e.assembly_id): e.id for e in session.query(Exon).all()}
    vdb_seq_syn_map = {(ss.sequence_region_id, ss.xref) for ss in session.query(SequenceRegionXref).all()}

    exon_xrefs = []

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if (vdb_exon_idx_map[idx], synonym, dataset_ids['CHESS']) not in vdb_seq_syn_map and synonym != '':
                    exon_xrefs.append(SequenceRegionXref(vdb_exon_idx_map[idx], synonym, 'synonym', dataset_ids['CHESS'], dataset_ids[dataset_name]))

    session.add_all(exon_xrefs)
    session.commit()


def load_chess_transcripts(session, transcript_exon_gff_df):
    """Load transcript data from CHESS database.

    This function processes and loads transcript annotations from CHESS,
    including transcript-exon relationships and transcript attributes.

    Args:
        session: SQLAlchemy session object
        transcript_exon_gff_df (pd.DataFrame): CHESS GFF dataframe with 
            transcript and exon information.
            Expected columns include:
            - ID: Transcript identifier
            - start_tx: Transcript start position
            - end_tx: Transcript end position
            - strand_tx: Strand orientation
            - assembly_id_tx: Assembly identifier
            - gene_name_tx: Associated gene name
            - gene_type_tx: Gene type
            - db_xref_tx: Database cross-references
            - num_samples_tx: Number of samples
            - max_tpm_tx: Maximum TPM value

    Returns:
        dict: Mapping of transcript IDs to their exon information
    """

    unique_cols = ['ID', 'start_tx', 'end_tx', 'strand_tx', 'assembly_id_tx', 'gene_name_tx', 'gene_type_tx', 'db_xref_tx', 'num_samples_tx', 'max_tpm_tx']
    transcript_exon_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_transcript_exon_gff_df = transcript_exon_gff_df.groupby(unique_cols).aggregate(list)

    vdb_transcript_idx_map = {t.transcript_idx: t.id for t in session.query(Transcript).all()}
    vdb_exon_map = {(e.start, e.end, e.strand, e.assembly_id): e.id for e in session.query(Exon).all()}

    transcripts = []
    transcript_exon_map = {}
    synonym_dict = {}

    existing_entries = []
    update_entries = []
    gene_update_entries = []

    unmapped_gene = session.query(Gene).filter(Gene.hgnc_id == 'unmapped').one()

    for (transcript_id, start, end, strand, assembly_id, gene_name, gene_type, db_xref, num_samples, max_tpm), row in grouped_transcript_exon_gff_df.iterrows():

        try:
            gene_id = synonym_vdb_gene_map[gene_name]
            gene_update_entries.append({"id": gene_id, "chess_id": transcript_id, 
                                        "attrs": {"gene_type": gene_type,
                                                  "db_xref": db_xref,
                                                  "num_samples": num_samples,
                                                  "max_tpm": max_tpm}})
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

            update_entries.append({"id": transcript_vdb_id, "chess_id": transcript_id, 
                                "attrs": {"gene_type": gene_type,
                                          "db_xref": db_xref,
                                          "num_samples": num_samples,
                                          "max_tpm": max_tpm}})

        elif transcript_idx not in existing_entries:

            existing_entries.append(transcript_idx)
            
            transcripts.append(Transcript(start, end, strand, assembly_id,
                                          block_sizes, chrom_starts, transcript_idx, transcript_idx_str,
                                          gene_id, '', '', transcript_id, '', 
                                          '', gene_type, 
                                          attrs={"gene_type": gene_type,
                                                 "db_xref": db_xref,
                                                 "num_samples": num_samples,
                                                 "max_tpm": max_tpm}))

        syn_dict = {source: set() for source in row.source_ex}

        for i, ID in enumerate(row.start_ex):
            syn_dict[row.source_ex[i]].add(row.db_xref_tx)
            syn_dict[row.source_ex[i]].add(transcript_id)

        synonym_dict[transcript_idx] = syn_dict

        transcript_exon_map[transcript_idx] = []

        for (exon_num, exon_id) in exon_ids:
            transcript_exon_map[transcript_idx].append((exon_num, exon_id))


    session.bulk_update_mappings(Gene, gene_update_entries,)
    session.commit()

    session.bulk_update_mappings(Transcript, update_entries,)
    session.commit()

    session.add_all(transcripts)
    session.commit()

    logging.info(f'Updated {len(gene_update_entries)} genes with CHESS info that had a synonym match to GENCODE/RefSeq')
    logging.info(f'Updated {len(update_entries)} transcripts with CHESS info that had exact coordinate matches to GENCODE/RefSeq')
    logging.info(f'Added {len(transcripts)} CHESS transcripts without an exact coordinate match to GENCODE/RefSeq')

    vdb_transcript_idx_map = {t.transcript_idx: t.id for t in session.query(Transcript).all()}

    transcript_xrefs = []

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                transcript_xrefs.append(TranscriptXref(vdb_transcript_idx_map[idx], synonym, 'synonym', 
                                                       dataset_ids['CHESS'], dataset_ids['CHESS']))

    session.add_all(transcript_xrefs)
    session.commit()

    return transcript_exon_map


def load_chess_transcript_exons(session, transcript_exon_map):
    """
    """
    vdb_transcript_map = {t.transcript_idx: t.id for t in session.query(Transcript).all()}
    
    transcript_exons = []

    existing_entries = set([(te.transcript_id, te.exon_id, te.exon_number) for te in session.query(TranscriptExon).all()])
    entries = set()

    for transcript_idx, vals in transcript_exon_map.items():
        for (exon_number, exon_id) in vals:
            entries.add((vdb_transcript_map[transcript_idx], exon_id, exon_number))

    entries.difference_update(existing_entries)
            
    for entry in entries:
        transcript_exons.append(TranscriptExon(*entry, attrs={'mapping_source': 'CHESS'}))
        
    session.add_all(transcript_exons)
    session.commit()

    logging.info(f'Added {len(transcript_exons)} transcript <-> exon mappings from CHESS')


def load_chess_cds(session, cds_gff_df):
    """Load CDS data from CHESS database.

    This function processes and loads coding sequence (CDS) annotations from CHESS,
    creating new CDS entries and their cross-references.

    Args:
        session: SQLAlchemy session object
        cds_gff_df (pd.DataFrame): CHESS GFF dataframe filtered for CDS regions.
            Expected columns include:
            - start: CDS start position
            - end: CDS end position
            - strand: Strand orientation
            - assembly_id: Assembly identifier
            - Parent: Parent transcript identifier
            - source: Data source
            - ID: CDS identifier
            - db_xref: Database cross-references

    Returns:
        None
    """

    unique_cols = ['start', 'end', 'strand', 'assembly_id']
    cds_gff_df.sort_values(by=unique_cols, inplace=True)
    grouped_cds_gff_df = cds_gff_df.groupby(unique_cols).aggregate(list)

    vdb_cds_idx_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}

    cds = []
    synonym_dict = {}
    update_entries = []
    existing_entries = []

    for cds_idx, row in grouped_cds_gff_df.iterrows():
        chess_id = f'cds-{row.Parent[0]}'

        if cds_idx in vdb_cds_idx_map.keys():
            cds_id = vdb_cds_idx_map[cds_idx]

            update_entries.append({"id": cds_id, "chess_id": chess_id})

        elif cds_idx not in existing_entries:

            existing_entries.append(cds_idx)
            
            cds.append(Cds(cds_idx[0], cds_idx[1], cds_idx[2], cds_idx[3],
                        '', '', chess_id, '', '', '', ''))

        syn_dict = {source: set() for source in row.source}

        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(f'cds-{row.Parent[i]}')
            syn_dict[row.source[i]].add(f'cds-{row.db_xref[i]}')

        synonym_dict[cds_idx] = syn_dict

    session.bulk_update_mappings(Cds, update_entries,)
    session.commit()

    session.add_all(cds)
    session.commit()

    logging.info(f'Updated {len(update_entries)} CDS with CHESS info that had exact coordinate matches to GENCODE/RefSeq')
    logging.info(f'Added {len(cds)} CHESS CDS without an exact coordinate match to GENCODE/RefSeq')

    vdb_cds_idx_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}
    vdb_seq_syn_map = {(ss.sequence_region_id, ss.xref) for ss in session.query(SequenceRegionXref).all()}

    cds_xrefs = []

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if (vdb_cds_idx_map[idx], synonym, dataset_ids['CHESS']) not in vdb_seq_syn_map and synonym != '':
                    cds_xrefs.append(SequenceRegionXref(vdb_cds_idx_map[idx], synonym, 'synonym', dataset_ids['CHESS'], dataset_ids[dataset_name]))

    session.add_all(cds_xrefs)
    session.commit()


def load_openprot_cds(session, bed_df):
    """Load CDS data from OpenProt database.

    This function processes and loads coding sequence (CDS) annotations from 
    OpenProt, handling both new entries and updates to existing ones.

    Args:
        session: SQLAlchemy session object
        bed_df (pd.DataFrame): OpenProt BED format dataframe.
            Expected columns include:
            - name: OpenProt identifier
            - chromStart: Start position
            - chromEnd: End position
            - strand: Strand orientation
            - assembly_id: Assembly identifier
            - blockCount: Number of blocks
            - blockSizes: Block sizes
            - blockStarts: Block start positions

    Returns:
        None

    Note:
        Only processes entries with names starting with 'IP' or 'II'
    """

    vdb_cds_idx_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}

    cds = []
    synonym_dict = {}
    update_entries = []
    existing_entries = set()

    for idx, row in bed_df.iterrows():
        
        if not row['name'][0:2] in ['IP', 'II']:
            continue
        
        block_sizes = [int(x) for x in row.blockSizes.split(',')]
        block_starts = [int(x) for x in row.blockStarts.split(',')]
        
        for i in range(row.blockCount):
            start = row.chromStart + block_starts[i] + 1
            end = row.chromStart + block_starts[i] + block_sizes[i]
            
            cds_idx = (start, end, row.strand, row.assembly_id)
            
            if cds_idx in vdb_cds_idx_map.keys():
                cds_id = vdb_cds_idx_map[cds_idx]

                update_entries.append({"id": cds_id, "openprot_id": f"cds-{row['name']}-{i}"})

            elif cds_idx not in existing_entries:

                existing_entries.add(cds_idx)
                cds.append(Cds(cds_idx[0], cds_idx[1], cds_idx[2], cds_idx[3],
                            '', '', '', '', '', '', '', f"cds-{row['name']}-{i}"))

            synonym_dict[cds_idx] = {"openprot": set([f"cds-{row['name']}-{i}"])}

    session.bulk_update_mappings(Cds, update_entries,)
    session.commit()

    session.add_all(cds)
    session.commit()

    logging.info(f'Updated {len(update_entries)} CDS with Openprot info that had exact coordinate matches to GENCODE/Refseq/CHESS')
    logging.info(f'Added {len(cds)} Openprot CDS without an exact coordinate match to GENCODE/Refseq/CHESS')

    vdb_cds_idx_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}

    cds_xrefs = []

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '': continue
                cds_xrefs.append(SequenceRegionXref(vdb_cds_idx_map[idx], synonym, 'synonym', dataset_ids[dataset_name], dataset_ids[dataset_name]))

    session.add_all(cds_xrefs)
    session.commit()


def load_openprot_orfs(session, bed_df):
    """Load ORF data from OpenProt database.

    This function processes and loads open reading frame (ORF) annotations from
    OpenProt, creating new entries and their relationships with CDS regions.

    Args:
        session: SQLAlchemy session object
        bed_df (pd.DataFrame): OpenProt BED format dataframe.
            Expected columns include:
            - name: OpenProt identifier
            - chromStart: Start position
            - chromEnd: End position
            - strand: Strand orientation
            - assembly_id: Assembly identifier
            - blockCount: Number of blocks
            - blockSizes: Block sizes
            - blockStarts: Block start positions
            - phases: Phase information for each block

    Returns:
        dict: Mapping of CDS IDs to their ORF information

    Note:
        Only processes entries with names starting with 'IP' or 'II'
    """
    vdb_cds_idx_map = {(c.start, c.end, c.strand, c.assembly_id): c.id for c in session.query(Cds).all()}
    vdb_orf_idx_map = {o.orf_idx: o.id for o in session.query(Orf).all()}

    orfs = []
    cds_orf_map = {}
    existing_entries = set()
    update_entries = []

    synonym_dict = {}

    for i, row in bed_df.iterrows():
        
        if not row['name'][0:2] in ['IP', 'II']:
            continue
            
        block_sizes = [int(x) - 1 for x in row.blockSizes.split(',')]
        block_starts = [int(x) for x in row.blockStarts.split(',')]
        chrom_starts = [row.chromStart + 1 + int(x) for x in row.blockStarts.split(',')]
        phases = ast.literal_eval(row.phases)
        cds_ids = []

        for i in range(row.blockCount):
            start = row.chromStart + block_starts[i] + 1
            end = row.chromStart + block_starts[i] + block_sizes[i] + 1
            cds_ids.append((i, phases[i], vdb_cds_idx_map[(start, end, row.strand, row.assembly_id)]))
        
        orf_start = row.chromStart + 1
        orf_end = row.chromEnd
        block_sizes = ';'.join([str(x) for x in block_sizes])
        chrom_starts = ';'.join([str(x) for x in chrom_starts])
        exon_frames = ';'.join([str(x) for x in phases])
        orf_val = f'{orf_start}_{orf_end}_{row.strand}_{row.assembly_id}_{block_sizes}_{chrom_starts}_{exon_frames}'
        orf_idx = hashlib.sha3_512(orf_val.encode('utf-8')).hexdigest()

        if orf_idx in vdb_orf_idx_map.keys():
            orf_id = vdb_orf_idx_map[orf_idx]

            update_entries.append({"id": orf_id, "openprot_id": row['name']})

        elif orf_idx not in existing_entries:

            existing_entries.add(orf_idx)

            orfs.append(Orf(orf_start, orf_end, row.strand, row.assembly_id, 
                        block_sizes, chrom_starts, exon_frames, orf_idx, orf_val,
                        '', '', '', '', '', '', {}, row['name']))
        
        synonym_dict[orf_idx] = {"openprot": set([row['name']])}

        
        cds_orf_map.update({cds_id: (orf_idx, exon_num, phase) \
                                        for (exon_num, phase, cds_id) in cds_ids})

    session.bulk_update_mappings(Orf, update_entries,)
    session.commit()

    session.add_all(orfs)
    session.commit()

    logging.info(f'Updated {len(update_entries)} ORFs with Openprot info that had exact coordinate matches to GENCODE/RefSeq/CHESS')
    logging.info(f'Added {len(orfs)} Openprot ORFs without an exact coordinate match to GENCODE/RefSeq/CHESS')

    orf_xrefs = []

    vdb_orf_idx_map = {o.orf_idx: o.id for o in session.query(Orf).all()}

    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '': continue
                orf_xrefs.append(OrfXref(vdb_orf_idx_map[idx], synonym, 'synonym', dataset_ids[dataset_name], dataset_ids[dataset_name]))

    session.add_all(orf_xrefs)
    session.commit()

    return cds_orf_map


def update_riboseq_cds(session, gencode_dir, version):
    """Update CDS information from Ribo-seq data.

    This function processes Ribo-seq CDS annotations, updating existing CDS entries 
    with corrected one-based coordinates and creating appropriate cross-references.
    It specifically handles coordinate conversion and synonym mapping for 
    ribosome profiling data.

    Args:
        session: SQLAlchemy session object
        gencode_dir (Path): Directory containing GENCODE files
        version (str): GENCODE version string (e.g., 'v42')

    Returns:
        None

    Note:
        The function performs two main tasks:
        1. Updates CDS coordinates to use one-based system
        2. Creates cross-references between CDS entries and their Ribo-seq identifiers

        The function expects a 'Ribo-seq_ORFs.bed' file in the GENCODE directory
        with the following columns:
        - chrom: Chromosome name
        - chromStart: Start position (0-based)
        - chromEnd: End position
        - name: Ribo-seq identifier
        - score: Score value
        - strand: Strand orientation ('+' or '-')
        - blockCount: Number of blocks
        - blockSizes: Comma-separated list of block sizes
        - chromStarts: Comma-separated list of block starts
        - phases: List of phase values
    """
    bed_cols = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 
                'strand', 'thickStart', 'thickEnd', 'itemRgb', 
                'blockCount', 'blockSizes', 'chromStarts', 'name2', 
                'cdsStartStat', 'cdsEndStat', 'exonFrames', 'type',
                'geneName', 'geneName2', 'geneType', 'transcript_biotype',
                'sequence', 'all_transcript_ids', 'all_gene_ids', 
                'replicated', 'ref_studies']

    orf_bed_df = pd.read_csv(
        gencode_dir.joinpath(version, 'Ribo-seq_ORFs.bed'),
        names=bed_cols,
        sep='\t'
    )

    vdb_cds_idx_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id 
        for c in session.query(Cds).all()
    }

    update_entries = []
    synonym_dict = {}

    for idx, row in orf_bed_df.iterrows():
        for i, (block, cstart) in enumerate(
            zip(row.blockSizes.rstrip(',').split(','),
                row.chromStarts.rstrip(',').split(',')), 
            start=1
        ):
            cds_start = int(row.chromStart) + int(cstart)
            cds_end = cds_start + int(block)

            cds_idx = (cds_start, cds_end, row.strand, assembly_ids[row.chrom])

            cds_id = vdb_cds_idx_map[cds_idx]
            
            ## Update coordinates to 1-based 
            chromStart = int(row.chromStart) + 1
            
            cds_start = int(chromStart) + int(cstart)
            cds_end = cds_start + int(block) - 1
            
            cds_idx = (cds_start, cds_end, row.strand, assembly_ids[row.chrom])
            
            if cds_idx in vdb_cds_idx_map:
                synonym_dict[cds_idx] = {"ENSEMBL": set([f"cds-{row['name']}-{i}"])}
            else: 
                update_entries.append({
                    "id": cds_id,
                    "start": cds_start,
                    "end": cds_end
                })

    session.bulk_update_mappings(Cds, update_entries)
    session.commit()

    logging.info(
        f'Updated {len(update_entries)} RiboSeq CDS with corrected '
        'one-based coordinates'
    )

    vdb_cds_idx_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id 
        for c in session.query(Cds).all()
    }

    cds_xrefs = []
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                if synonym == '':
                    continue
                cds_xrefs.append(SequenceRegionXref(
                    vdb_cds_idx_map[idx],
                    synonym,
                    'synonym',
                    dataset_ids[dataset_name],
                    dataset_ids[dataset_name]
                ))

    session.add_all(cds_xrefs)
    session.commit()

    logging.info(
        f'Added {len(cds_xrefs)} RiboSeq CDS synonyms after coordinate '
        'update to one-based'
    )


def load_psl_phase_cds(session, psl_df, source):
    """Load CDS information from PSL data.

    Args:
        session: SQLAlchemy session object
        psl_df (pd.DataFrame): PSL alignment dataframe
        source (str): Data source identifier
    """
    cds = []
    synonym_dict = {}

    for i, row in psl_df.iterrows():
        cds.append(Cds(
            row.tStart,
            row.tEnd,
            row.strand,
            row.assembly_id,
            row.qName,
            '', '', '',
            '', '', source
        ))

        synonym_dict[(
            row.tStart,
            row.tEnd,
            row.strand,
            row.assembly_id
        )] = {'PSL': set([row.qName])}

    session.add_all(cds)
    session.commit()

    logging.info(f'Added {len(cds)} {source} CDS')

    cds_xrefs = []
    ensembl_orf_cds_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id
        for c in session.query(Cds).all()
    }
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                cds_xrefs.append(SequenceRegionXref(
                    ensembl_orf_cds_map[idx],
                    synonym,
                    'synonym',
                    dataset_ids['PSL'],
                    dataset_ids[dataset_name]
                ))

    session.add_all(cds_xrefs)
    session.commit()

    logging.info(f'Added {len(cds_xrefs)} {source} CDS synonyms')


def load_psl_phase_orfs_legacy(session, psl_df, phase_df, source):
    """Load legacy ORF information from PSL and phase data.

    Args:
        session: SQLAlchemy session object
        psl_df (pd.DataFrame): PSL alignment dataframe
        phase_df (pd.DataFrame): Phase information dataframe
        source (str): Data source identifier
    """
    orfs = []
    cds_orf_map = {}
    synonym_dict = {}
    ensembl_orf_cds_map = {
        (c.start, c.end, c.strand, c.assembly_id): c.id
        for c in session.query(Cds).all()
    }

    for i, row in psl_df.iterrows():
        block_sizes = [str(x) for x in row.blockSizes]
        chrom_starts = [str(x) for x in row.tStarts]
        exon_frames = [str(x) for x in row.qStarts]

        block_sizes = ';'.join(block_sizes)
        chrom_starts = ';'.join(chrom_starts)
        exon_frames = ';'.join(exon_frames)

        orf_val = (
            f'{row.tStart}_{row.tEnd}_{row.strand}_{row.assembly_id}_'
            f'{block_sizes}_{chrom_starts}_{exon_frames}'
        )
        orf_idx = hashlib.sha3_512(orf_val.encode('utf-8')).hexdigest()

        if orf_idx in synonym_dict.keys():
            synonym_dict[orf_idx]['PSL'].add(row.qName)
            continue
        else:
            synonym_dict[orf_idx] = {'PSL': set([row.qName])}

        orfs.append(Orf(
            row.tStart,
            row.tEnd,
            row.strand,
            row.assembly_id,
            block_sizes,
            chrom_starts,
            exon_frames,
            orf_idx,
            orf_val,
            '', '', '',
            row.qName,
            '', source
        ))

        cds_orf_map[orf_idx] = []
        cds_orf_map[orf_idx].append((
            1, 0,
            ensembl_orf_cds_map[
                (row.tStart,
                 row.tEnd,
                 row.strand,
                 row.assembly_id)
            ]
        ))

    session.add_all(orfs)
    session.commit()

    logging.info(f'Added {len(orfs)} {source} ORFs')

    orf_xrefs = []
    ensembl_orf_orf_map = {
        o.orf_idx: o.id
        for o in session.query(Orf).all()
    }
    dataset_ids = {x.name: x.id for x in session.query(Dataset).all()}

    for idx, synonym_entries in synonym_dict.items():
        for dataset_name, synonyms in synonym_entries.items():
            for synonym in synonyms:
                orf_xrefs.append(SequenceRegionXref(
                    ensembl_orf_orf_map[idx],
                    synonym,
                    'synonym',
                    dataset_ids['PSL'],
                    dataset_ids[dataset_name]
                ))

    session.add_all(orf_xrefs)
    session.commit()

    logging.info(f'Added {len(orf_xrefs)} {source} ORF synonyms')

    return cds_orf_map


def load_psl_phase_proteins(session, phase_df, source):
    """Load protein information from phase data.

    Args:
        session: SQLAlchemy session object
        phase_df (pd.DataFrame): Phase information dataframe
        source (str): Data source identifier
    """
    proteins = []
    for i, row in phase_df.iterrows():
        proteins.append(Protein(
            row.id,
            row.id,
            'active',
            row.id,
            '',
            'Homo sapiens',
            len(row.aa),
            row.aa,
            '', '', '', '', '', '', '', '', '', '', '', '', ''
        ))

    session.add_all(proteins)
    session.commit()

    logging.info(f'Added {len(proteins)} {source} proteins')


def update_psl_phase_orfs(session, psl_df):
    """Update ORF information from PSL data.

    Args:
        session: SQLAlchemy session object
        psl_df (pd.DataFrame): PSL alignment dataframe
    """
    update_entries = []
    ensembl_orf_orf_map = {
        o.ensembl_id: o.id
        for o in session.query(Orf).all()
    }

    for i, row in psl_df.iterrows():
        if row.qName not in ensembl_orf_orf_map.keys():
            continue

        block_sizes = [str(x) for x in row.blockSizes]
        chrom_starts = [str(x) for x in row.tStarts]
        exon_frames = [str(x) for x in row.qStarts]

        block_sizes = ';'.join(block_sizes)
        chrom_starts = ';'.join(chrom_starts)
        exon_frames = ';'.join(exon_frames)

        orf_val = (
            f'{row.tStart}_{row.tEnd}_{row.strand}_{row.assembly_id}_'
            f'{block_sizes}_{chrom_starts}_{exon_frames}'
        )
        orf_idx = hashlib.sha3_512(orf_val.encode('utf-8')).hexdigest()

        update_entries.append({
            'id': ensembl_orf_orf_map[row.qName],
            'start': row.tStart,
            'end': row.tEnd,
            'strand': row.strand,
            'assembly_id': row.assembly_id,
            'block_sizes': block_sizes,
            'chrom_starts': chrom_starts,
            'exon_frames': exon_frames,
            'orf_idx': orf_idx,
            'orf_val': orf_val
        })

    session.bulk_update_mappings(Orf, update_entries)
    session.commit()

    logging.info(f'Updated {len(update_entries)} ORFs')


def load_gencode_lncRNA_genes(session, gene_gff_df):
    """Load long non-coding RNA genes from GENCODE GFF data.

    This function processes GENCODE gene annotations specifically for lncRNA genes,
    adding them to the database and creating appropriate cross-references.
    It handles both new gene entries and updates to existing genes.

    Args:
        session: SQLAlchemy session object
        gene_gff_df (pd.DataFrame): GENCODE GFF dataframe filtered for lncRNA genes.
            Expected columns include:
            - start: Gene start position
            - end: Gene end position
            - strand: Strand orientation ('+' or '-')
            - assembly_id: Reference assembly identifier
            - ID: GENCODE gene identifier
            - gene_type: Type of gene (e.g., 'lncRNA')
            - tag: Additional gene annotations
            - gene_name: Official gene symbol
            - hgnc_id: HGNC identifier
            - source: Data source identifier

    Returns:
        None

    Note:
        The function creates both gene entries and cross-references (synonyms)
        in the database. It handles duplicate entries by updating existing records
        rather than creating new ones.
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
        if gene_idx in vdb_gene_idx_map.keys():
            gene_id = vdb_gene_idx_map[gene_idx]
            update_entries.append({
                "id": gene_id,
                "ensembl_id": row.ID[0],
                "attrs": {
                    "gene_type": row.gene_type[0],
                    "tag": row.tag[0]
                }
            })
        elif gene_idx not in synonym_dict.keys():
            new_genes.append(Gene(
                gene_idx[0], gene_idx[1], gene_idx[2], gene_idx[3],
                row.hgnc_id[0], row.gene_name[0], row.ID[0], '', '', '',
                row.gene_type[0], '',
                attrs={"gene_type": row.gene_type[0], "tag": row.tag[0]}
            ))

        syn_dict = {source: set() for source in row.source}
        for i, ID in enumerate(row.ID):
            syn_dict[row.source[i]].add(row.gene_name[i])
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
                if synonym == '':
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