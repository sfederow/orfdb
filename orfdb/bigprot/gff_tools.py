"""
This module provides tools for working with GFF (General Feature Format) files. It includes functions to download GFF files, load them into pandas DataFrames, 
filter and find attributes within these DataFrames, and sparsify the data for easier handling. Additionally, it contains functions for working with GTF (Gene Transfer Format) files, 
including parsing attributes and sparsifying the data. The module also includes functionality to generate transcripts from GTF data.
"""
import os
import json
import subprocess
from tqdm import tqdm
from typing import Dict, List, Set, Any

import numpy as np
import pandas as pd

from . import (constants,
               misc,
               coordinates)


def download_gff(output_folder: str, source_url: str) -> None:
    """
    Downloads a GFF file from a specified URL into an output folder.

    Args:
        output_folder (str): The folder where the GFF file will be saved.
        source_url (str): The URL from which the GFF file will be downloaded.
    """
    os.makedirs(output_folder, exist_ok=True)
    cmd = f'wget -nd -P {output_folder} -r -l 1 {source_url}'
    print(subprocess.check_output(cmd, shell=True,
          stderr=subprocess.STDOUT, encoding='utf-8'))


def load_gff(gff_fpath: str) -> pd.DataFrame:
    """
    Loads a GFF file into a pandas DataFrame.

    Args:
        gff_fpath (str): The file path to the GFF file.

    Returns:
        pd.DataFrame: The loaded GFF data as a pandas DataFrame.
    """
    gff_df = pd.read_csv(gff_fpath, sep='\t', comment='#',
                         header=None, names=constants.GFF_FIELDS, compression='gzip')

    parsed_attributes = []
    for entry in tqdm(gff_df.attributes):
        splat = entry.split(';')

        attrs = dict([element.split('=') for element in splat])
        parsed_attributes.append(attrs)

    gff_df['attributes'] = parsed_attributes
    return gff_df


def filter_by(df: pd.DataFrame, col: str, include_values: List[str] = None, exclude_values: List[str] = None) -> pd.DataFrame:
    """
    Filters a DataFrame by including or excluding values in a specified column.

    Args:
        df (pd.DataFrame): The DataFrame to filter.
        col (str): The column name to filter by.
        include_values (List[str], optional): The values to include. Defaults to None.
        exclude_values (List[str], optional): The values to exclude. Defaults to None.

    Returns:
        pd.DataFrame: The filtered DataFrame.
    """
    if include_values is None:
        include_values = []
    if exclude_values is None:
        exclude_values = []
    return df.loc[[val in set(include_values) and val not in set(exclude_values) for val in df[col]]]


def find_all_attributes(gff_df: pd.DataFrame, attribute_column: str = 'attribute') -> Set[str]:
    """
    Finds all unique attributes in a specified column of a GFF DataFrame.

    Args:
        gff_df (pd.DataFrame): The GFF DataFrame.
        attribute_column (str, optional): The column name to search for attributes. Defaults to 'attribute'.

    Returns:
        Set[str]: A set of all unique attributes found.
    """
    all_attrs = set([])
    for row_idx in tqdm(gff_df.index):
        all_attrs.update(gff_df.loc[row_idx, attribute_column].keys())
    return all_attrs


def sparsify_dicts(dense_df: pd.DataFrame, columns: List[str], null_value: str = '') -> pd.DataFrame:
    """
    Converts dense dictionary columns in a DataFrame to sparse format.

    Args:
        dense_df (pd.DataFrame): The dense DataFrame.
        columns (List[str]): The columns to sparsify.
        null_value (str, optional): The value to use for missing data. Defaults to ''.

    Returns:
        pd.DataFrame: The DataFrame with sparsified columns.
    """
    new_cols = set([])
    for col in columns:
        new_cols.update(find_all_attributes(dense_df, col))

    new_columns = {col: {} for col in new_cols}

    for col in columns:
        for row_idx, attrs in tqdm(dense_df[col].items()):
            for k, v in attrs.items():
                new_columns[k][row_idx] = v

    sparse_df = dense_df.copy()
    for col, col_values in new_columns.items():
        sparse_df[col] = pd.Series(col_values)
        sparse_df[col] = sparse_df[col].fillna(null_value)

    return sparse_df


def find_all_tags(gff_df: pd.DataFrame, tag_column: str, sep: str = ',') -> Set[str]:
    """
    Finds all unique tags in a specified column of a DataFrame, where tags are separated by a specified separator.

    Args:
        gff_df (pd.DataFrame): The DataFrame to search.
        tag_column (str): The column name to search for tags.
        sep (str, optional): The separator used between tags. Defaults to ','.

    Returns:
        Set[str]: A set of all unique tags found.
    """
    all_tags = set([])
    for tag_entry in tqdm(gff_df[tag_column]):
        all_tags.update(tag_entry.split(sep))
    return all_tags


def sparsify_lists(dense_df: pd.DataFrame, columns: List[str], sep: str = ',') -> pd.DataFrame:
    """
    Converts dense list columns in a DataFrame to sparse format, where lists are separated by a specified separator.

    Args:
        dense_df (pd.DataFrame): The dense DataFrame.
        columns (List[str]): The columns to sparsify.
        sep (str, optional): The separator used between list items. Defaults to ','.

    Returns:
        pd.DataFrame: The DataFrame with sparsified columns.
    """
    new_cols = set([])

    for col in columns:
        new_cols.update(find_all_tags(dense_df, col, sep=sep))

    new_columns = {col: {} for col in new_cols}

    for col in columns:
        for row_idx, cell_value in tqdm(zip(dense_df.index, dense_df[col])):
            tags = cell_value.split(sep)
            for tag in tags:
                new_columns[tag][row_idx] = True

    col_series = {}
    for col, col_values in new_columns.items():
        col_series[col] = pd.Series(col_values, dtype=bool)

    null_value = False
    new_col_df = pd.DataFrame(col_series).fillna(null_value)

    return pd.concat((dense_df, new_col_df), axis=1)


def split_gtf_attributes(attribute_string: str, sep: str = ';') -> Dict[str, Any]:
    """
    Splits a GTF attribute string into a dictionary.

    Args:
        attribute_string (str): The attribute string to split.
        sep (str, optional): The separator used between attributes. Defaults to ';'.

    Returns:
        Dict[str, Any]: A dictionary of attributes.
    """
    attributes_dict = {}
    for element in attribute_string.strip().rstrip(sep).split(sep):
        element = element.strip()
        if '{' in element and '}' in element:
            assert 'dict' not in attributes_dict, attribute_string
            this_dict = json.loads(element)
            attributes_dict['dict'] = this_dict
            print(attributes_dict)
            print()
        hinge = element.find(' ')
        key, value = element[:hinge], element[hinge+1:].strip('"')
        attributes_dict[key] = value
    return attributes_dict


def find_all_attributes_gtf(gtf_df: pd.DataFrame, attribute_column: str = 'attribute') -> Set[str]:
    """
    Finds all unique attributes in a specified column of a GTF DataFrame.

    Args:
        gtf_df (pd.DataFrame): The GTF DataFrame.
        attribute_column (str, optional): The column name to search for attributes. Defaults to 'attribute'.

    Returns:
        Set[str]: A set of all unique attributes found.
    """
    all_attrs = set([])
    for attribute_string in tqdm(gtf_df[attribute_column]):
        all_attrs.update(split_gtf_attributes(attribute_string).keys())
    return all_attrs


def sparsify_gtf(dense_df: pd.DataFrame, columns: List[str] = ['attribute'], null_value: str = '', sep: str = ';') -> pd.DataFrame:
    """
    Converts dense GTF columns in a DataFrame to sparse format.

    Args:
        dense_df (pd.DataFrame): The dense DataFrame.
        columns (List[str], optional): The columns to sparsify. Defaults to ['attribute'].
        null_value (str, optional): The value to use for missing data. Defaults to ''.
        sep (str, optional): The separator used between attributes. Defaults to ';'.

    Returns:
        pd.DataFrame: The DataFrame with sparsified columns.
    """
    new_cols = set([])
    for col in columns:
        new_cols.update(find_all_attributes_gtf(dense_df, col))

    new_columns = {col: {} for col in new_cols}

    for col in columns:
        for row_idx, attribute_string in tqdm(zip(dense_df.index, dense_df[col])):
            for k, v in split_gtf_attributes(attribute_string).items():
                new_columns[k][row_idx] = v

    sparse_df = dense_df.copy()
    for col, col_values in new_columns.items():
        sparse_df[col] = pd.Series(col_values)
        sparse_df[col] = sparse_df[col].fillna(null_value)

    return sparse_df


def load_gtf(gtf_fpath: str) -> pd.DataFrame:
    return pd.read_csv(gtf_fpath, sep='\t', names=constants.GFF_FIELDS,
                       low_memory=False, comment='#').sort_values(['seqname', 'start'])


def load_and_preprocess_gtf(gtf_fpath: str) -> pd.DataFrame:
    """
    Loads a GTF file and preprocesses it by filtering for exons and sparsifying the data.

    Args:
        gtf_fpath (str): The file path to the GTF file.

    Returns:
        pd.DataFrame: The preprocessed GTF data as a pandas DataFrame.
    """

    sparse_exon_gtf = sparsify_gtf(
        filter_by(load_gtf(gtf_fpath=gtf_fpath), 'feature', ['exon']))

    return sparse_exon_gtf


def generate_transcripts_from_gtf(sparse_exon_gtf: pd.DataFrame, genome: Dict[str, str], accession_namespace=constants.DEFAULT_ACCESSION_NAMESPACE) -> List[Dict[str, Any]]:
    """
    Generates transcripts from a sparsified GTF DataFrame.

    Args:
        sparse_exon_gtf (pd.DataFrame): The sparsified GTF DataFrame containing exon data.
        genome (Dict[str, str]): A dictionary of genome sequences keyed by sequence name.

    Returns:
        List[Dict[str, Any]]: A list of dictionaries, each representing a transcript.
    """
    # ToDo: include transcripts in filtered gtf and validate them against the ones we reconstruct from exons.
    seq_lengths = {seqname: len(seq) for seqname, seq in genome.items()}

    transcripts = []

    for transcript_id, exons in sparse_exon_gtf.groupby('transcript_id'):
        chrom_starts = []
        chrom_ends = []
        block_sizes = []
        transcript_strand = None

        for exon_seqname, exon_strand, exon_start, exon_end in zip(exons.seqname, exons.strand, exons.start, exons.end):
            if transcript_strand is None:
                transcript_strand = exon_strand
            else:
                assert transcript_strand == exon_strand

            this_exon_interval = coordinates.GffInterval(
                start=exon_start, end=exon_end)
            chrom_starts.append(this_exon_interval.gff_start)
            block_sizes.append(this_exon_interval.gff_size)
            chrom_ends.append(this_exon_interval.gff_end)

        transcript_start, transcript_end = chrom_starts[0], chrom_ends[-1]

        this_transcript = {'assembly.id': exon_seqname,
                           f'assembly.{accession_namespace}_accession': exon_seqname,
                           'assembly.sequence_length': seq_lengths[exon_seqname],
                           'transcript.block_sizes': misc.encode_numeric_vector(block_sizes),
                           'transcript.chrom_starts': misc.encode_numeric_vector(chrom_starts),
                           'transcript.end': transcript_end,
                           'transcript.ensembl_id': transcript_id,
                           'transcript.id': transcript_id,
                           'transcript.start': transcript_start,
                           'transcript.strand': transcript_strand,
                           'transcript.transcript_type': 'protein_coding'}
        transcripts.append(this_transcript)
    return transcripts
