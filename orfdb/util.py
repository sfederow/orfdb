"""Utility functions for the ORF database."""

from orfdb.base import *
from orfdb import settings
from conservation import phylocsf

import ast
import re
import os
import logging
import pandas as pd
from pathlib import Path
from time import time
from sys import stdout


def get_or_create(session, query_class, **kwargs):
    """Query the database for an existing entry or create a new one.

    Args:
        session: SQLAlchemy session object
        query_class: The SQLAlchemy model class to query
        **kwargs: Keyword arguments used for filtering the query

    Returns:
        tuple: (query_result, bool_exists) where bool_exists is True if the
               entry already existed
    """
    res = session.query(query_class).filter_by(**kwargs).first()
    if res is not None:
        return res, True
    res = query_class(**kwargs)
    session.add(res)
    session.commit()
    return res, False


def timing(function):
    """Decorator to log execution time of functions.

    Args:
        function: The function to be timed

    Returns:
        wrapper: The wrapped function that includes timing functionality
    """
    def wrapper(*args, **kwargs):
        arg_str = str(args)
        if arg_str[-2] == ',':  # trailing comma
            arg_str = arg_str[:-2] + ')'
        try:
            name = function.__name__
        except AttributeError:
            name = function.func_name
        logging.debug('starting %s' % name)
        stdout.flush()
        start = time()
        res = function(*args, **kwargs)
        logging.debug('%s complete (%.2f sec)' % (name, time() - start))
        return res
    return wrapper


def get_entrez_id(row):
    """Extract Entrez Gene ID from a row's Dbxref field.

    Args:
        row: DataFrame row containing a Dbxref field

    Returns:
        str: The Entrez Gene ID if found, empty string otherwise
    """
    entrez_gene_id = ''
    for val in row.Dbxref.split(','):
        vals = val.split(':')
        if vals[0] == 'GeneID':
            entrez_gene_id = vals[1]
    return entrez_gene_id


def get_hgnc_id(row):
    """Extract HGNC ID from a row's Dbxref field.

    Args:
        row: DataFrame row containing a Dbxref field

    Returns:
        str: The HGNC ID if found, empty string otherwise
    """
    hgnc_id = ''
    if re.match('.*HGNC:HGNC:(\d+).*', row.Dbxref):
        hgnc_id = re.match('.*HGNC:HGNC:(\d+).*', row.Dbxref).groups()[0]
    return hgnc_id


def calc_cds_phase_frame_bed(row, chrom_length):
    """Calculate CDS phase and reading frame for BED format entries.

    Args:
        row: DataFrame row containing BED format fields
        chrom_length: Length of the chromosome

    Returns:
        tuple: (phases, reading_frames) Lists containing phase and frame values
               for each block
    """
    block_sizes = [int(b) for b in row.blockSizes.rstrip(',').split(',')]
    block_starts = [int(b) for b in row.blockStarts.rstrip(',').split(',')]
    chrom_starts = [row.chromStart + int(bs) for bs in block_starts]
    block_coords = [(row.chromStart + int(bs),
                    row.chromStart + int(bs) + int(block_sizes[i]))
                   for i, bs in enumerate(block_starts)]

    reading_frames = []
    phases = []
    intron_lengths = [block_coords[i+1][0] - block_coords[i][1]
                     for i in range(len(block_coords)-1)]

    if row.strand == '+':
        for i, block in enumerate(block_coords):
            if i == 0:
                phases.append(0)
                reading_frames.append(block[0] % 3)
            else:
                reading_frames.append(
                    (reading_frames[i-1] + intron_lengths[i-1]) % 3)
                phases.append(
                    (3 - ((block_sizes[i-1] - phases[i-1]) % 3)) % 3)

    elif row.strand == '-':
        block_coords.reverse()
        block_sizes.reverse()
        intron_lengths = [block_coords[i][0] - block_coords[i+1][1]
                         for i in range(len(block_coords)-1)]

        for i, block in enumerate(block_coords):
            if i == 0:
                phases.append(0)
                reading_frames.append((chrom_length - block[1]) % 3)
            else:
                reading_frames.append(
                    (reading_frames[i-1] + intron_lengths[i-1]) % 3)
                phases.append(
                    (3 - ((block_sizes[i-1] - phases[i-1]) % 3)) % 3)

        reading_frames.reverse()
        phases.reverse()

    return phases, reading_frames


def calc_cds_phase_frame_psl(row):
    """Calculate CDS phase and reading frame for PSL format entries.

    Args:
        row: DataFrame row containing PSL format fields

    Returns:
        tuple: (phases, reading_frames) Lists containing phase and frame values
               for each block
    """
    block_sizes = [int(b) for b in row.blockSizes.rstrip(',').split(',')]
    chrom_starts = [int(ts) for ts in row.tStarts.rstrip(',').split(',')]
    block_coords = [(chrom_starts[i], chrom_starts[i] + int(bs))
                   for i, bs in enumerate(block_sizes)]

    reading_frames = []
    phases = []

    if row.strand == '+':
        intron_lengths = [block_coords[i+1][0] - block_coords[i][1]
                         for i in range(len(block_coords)-1)]

        for i, block in enumerate(block_coords):
            if i == 0:
                phases.append(0)
                reading_frames.append(block[0] % 3)
            else:
                reading_frames.append(
                    (reading_frames[i-1] + intron_lengths[i-1]) % 3)
                phases.append(
                    (3 - ((block_sizes[i-1] - phases[i-1]) % 3)) % 3)

    elif row.strand == '-':
        block_coords.reverse()
        block_sizes.reverse()
        intron_lengths = [block_coords[i][0] - block_coords[i+1][1]
                         for i in range(len(block_coords)-1)]

        for i, block in enumerate(block_coords):
            if i == 0:
                phases.append(0)
                reading_frames.append((row.tSize - block[1]) % 3)
            else:
                reading_frames.append(
                    (reading_frames[i-1] + intron_lengths[i-1]) % 3)
                phases.append(
                    (3 - ((block_sizes[i-1] - phases[i-1]) % 3)) % 3)

        reading_frames.reverse()
        phases.reverse()

    return phases, reading_frames


def combine_dicts(dict1, dict2):
    """Combine two dictionaries with set values.

    Args:
        dict1: First dictionary with set values
        dict2: Second dictionary with set values

    Returns:
        dict: Combined dictionary where values for common keys are unions of sets
    """
    common_keys = set(dict1.keys()) & set(dict2.keys())

    combined_dict = {}
    for key in common_keys:
        combined_dict[key] = dict1[key].union(dict2[key])

    for key in set(dict1.keys()) - common_keys:
        combined_dict[key] = dict1[key]

    for key in set(dict2.keys()) - common_keys:
        combined_dict[key] = dict2[key]

    return combined_dict


def parse_array(s):
    """Parse string representation of arrays into Python lists or numpy arrays.

    Args:
        s: String representation of an array

    Returns:
        list/array: Parsed array or original string if parsing fails
    """
    try:
        return ast.literal_eval(s)
    except (ValueError, SyntaxError):
        return s


def etl_phase1to5_psl(session, velia_dir):
    """Extract, transform, and load phase 1-5 PSL data.

    Args:
        session: SQLAlchemy session object
        velia_dir: Directory containing phase 1-5 data files

    Returns:
        tuple: (psl_df, phase_df) DataFrames containing processed PSL and
               phase data
    """
    phase_df = pd.read_csv(
        velia_dir.joinpath('phase1to5_ids_harmonized_20230314.csv'))
    psl_file = velia_dir.joinpath('phase1to5_ids_harmonized_nt.pslx')
    psl_raw_df = phylocsf.load_blat_psl(psl_file, 'nt')

    assembly_ids = {}
    for assembly in session.query(Assembly).all():
        if assembly.ucsc_style_name.startswith('chr') and len(
                assembly.ucsc_style_name) < 6:
            assembly_ids[assembly.ucsc_style_name] = assembly
        else:
            assembly_ids[assembly.genbank_accession] = assembly

    psl_df = psl_raw_df.merge(phase_df, left_on='qName', right_on='id')
    psl_df['assembly_id'] = psl_df.apply(
        lambda x: assembly_ids[x.tName].id, axis=1)

    return psl_df, phase_df


def etl_phase6_psl(session, velia_dir):
    """Extract, transform, and load phase 6 PSL data.

    Args:
        session: SQLAlchemy session object
        velia_dir: Directory containing phase 6 data files

    Returns:
        tuple: (psl_df, phase_df) DataFrames containing processed PSL and
               phase data
    """
    phase_6_dir = velia_dir.joinpath('phase_6')

    phase_df = pd.read_csv(phase_6_dir.joinpath('phase6_construct_nt_seq.csv'))
    phase_df.rename(columns={
        'primary_id': 'id',
        'source': 'source_long',
        'endogenous_aa_seq': 'aa',
        'endogenous_nt_seq': 'nucl'
    }, inplace=True)
    phase_df['ProteinID'] = phase_df.apply(
        lambda x: f"sORF{x.id.split('_')[1]}{x.id.split('_')[2]}", axis=1)

    genscript_df = pd.read_excel(phase_6_dir.joinpath('phase6_standardized.xlsx'))
    genscript_df['id'] = genscript_df.apply(
        lambda x: '_'.join(x['name'].split('_')[0:3]), axis=1)
    genscript_df.rename(columns={'Order ID': 'Genscript ID'}, inplace=True)
    genscript_df = genscript_df.groupby('id').first()
    genscript_df.reset_index(inplace=True)

    phase_df = phase_df.merge(genscript_df, left_on='id', right_on='id')

    source_map = {
        'EMBL Catalog GWAS Autoimmune SNPs in CDS of OpenProt sORFs':
            'autoimmune_gwas',
        'In-house Plasma MassSpec': 'plasma_mass_spec',
        'Viral sORF literature search': 'viral_sORF',
        'mappable;Public MassSpec sORF search': 'public_mass_spec'
    }

    phase_df['source'] = phase_df.apply(
        lambda x: source_map[x.source_long], axis=1)

    for i, row in phase_df[['source', 'source_long']].drop_duplicates().iterrows():
        dataset_name = f"velia_phase6_{row.source}"
        session.upsert(
            Dataset,
            name=dataset_name,
            description=row.source_long,
            type="dataset",
            attrs={}
        )

    non_viral_source = [
        'EMBL Catalog GWAS Autoimmune SNPs in CDS of OpenProt sORFs',
        'In-house Plasma MassSpec',
        'mappable;Public MassSpec sORF search'
    ]

    non_viral_df = phase_df[phase_df['source_long'].isin(non_viral_source)]

    assembly_ids = {}
    for assembly in session.query(Assembly).all():
        if assembly.ucsc_style_name.startswith('chr') and len(
                assembly.ucsc_style_name) < 6:
            assembly_ids[assembly.ucsc_style_name] = assembly
        else:
            assembly_ids[assembly.genbank_accession] = assembly

    dataset = 'phase6_construct_nt_seq_non-viral'
    file_type = 'pslx'

    psl_file = Path(f'{phase_6_dir}/{dataset}.{file_type}')
    psl_raw_df = phylocsf.load_blat_psl(psl_file, 'nt')

    psl_df = psl_raw_df.merge(non_viral_df, left_on='qName', right_on='id')
    psl_df['assembly_id'] = psl_df.apply(
        lambda x: assembly_ids[x.tName].id, axis=1)

    psl_df['Phase'] = 6
    phase_df['Phase'] = 6

    return psl_df, phase_df