"""Utility functions for the ORF database."""

from orfdb.base import *
from orfdb import settings

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
    if re.match(r".*HGNC:HGNC:(\d+).*", row.Dbxref):
        hgnc_id = re.match(r".*HGNC:HGNC:(\d+).*", row.Dbxref).groups()[0]
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

