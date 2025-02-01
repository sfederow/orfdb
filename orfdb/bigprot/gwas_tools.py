"""
This module contains functions for loading, filtering, and processing SNP associations. 
It includes functions to load SNP associations from a file, filter SNP associations based on valid chromosomes, 
generate interval trees for SNP associations, and find SNPs near a given transcript.
"""
import logging
from typing import List, Dict, Union

import intervaltree
import numpy as np
import pandas as pd
from intervaltree import IntervalTree

from . import constants
from . import validation

logger = logging.getLogger(__name__)


def load_snp_associations(ebl_associations_fpath: str) -> pd.DataFrame:
    """
    Loads SNP associations from a file.

    Args:
        ebl_associations_fpath (str): The file path of the SNP associations.

    Returns:
        pd.DataFrame: A DataFrame containing the SNP associations.
    """
    logger.info('Loading SNP associations from %s ...' %
                ebl_associations_fpath)
    snp_df = pd.read_csv(ebl_associations_fpath, sep='\t', low_memory=False)
    initial_size = snp_df.shape[0]

    logger.info('Loaded %s SNP associations.' % initial_size)
    return snp_df


def filter_snp_associations(snp_df: pd.DataFrame,
                            valid_chromosomes: List[str] = constants.DEFAULT_VALID_CHROMOSOMES) -> pd.DataFrame:
    """
    Filters SNP associations based on valid chromosomes.

    Args:
        snp_df (pd.DataFrame): The DataFrame containing the SNP associations.
        valid_chromosomes (List[str], optional): A list of valid chromosomes. Defaults to constants.DEFAULT_VALID_CHROMOSOMES.

    Returns:
        pd.DataFrame: A DataFrame containing the filtered SNP associations.
    """
    logger.info('Filtering SNP associations based on valid chromosomes...')
    snp_df = snp_df.loc[np.isin(snp_df.CHR_ID, valid_chromosomes)]
    filtered_size = snp_df.shape[0]
    logger.info('%s SNP associations remain after removing invalid chromosome entries' %
                filtered_size)

    return snp_df


def generate_snp_interval_trees(snp_df: pd.DataFrame) -> Dict[str, intervaltree.IntervalTree]:
    """
    Generates interval trees for SNP associations.

    Args:
        snp_df (pd.DataFrame): The DataFrame containing the SNP associations.

    Returns:
        Dict[str, intervaltree.IntervalTree]: A dictionary of interval trees for SNP associations.
    """
    logger.info('Generating SNP interval trees ...')
    snp_intervals: dict[str, IntervalTree] = {}
    for chrom, chrom_snps in snp_df.groupby('CHR_ID'):
        snp_intervals[str(chrom)] = intervaltree.IntervalTree.from_tuples(zip(chrom_snps.CHR_POS.astype(int),
                                                                              chrom_snps.CHR_POS.astype(
                                                                                  int) + 1,
                                                                              chrom_snps.SNP_ID_CURRENT.astype(str)))
    return snp_intervals


def find_snps_near_transcript(transcript_dict: Dict[str, Union[str, int]],
                              snp_intervals: Dict[str, intervaltree.IntervalTree],
                              promoter_upstream_window_size: int = constants.DEFAULT_PROMOTER_UPSTREAM_WINDOW_SIZE) -> \
        List[str]:
    """
    Finds SNPs near a given transcript.

    Args:
        transcript_dict (Dict[str, Union[str, int]]): The dictionary containing the transcript information.
        snp_intervals (Dict[str, intervaltree.IntervalTree]): The dictionary of interval trees for SNP associations.
        promoter_upstream_window_size (int, optional): The size of the promoter upstream window. Defaults to constants.DEFAULT_PROMOTER_UPSTREAM_WINDOW_SIZE.

    Returns:
        List[str]: A list of rs IDs for the SNPs found near the transcript.
    """
    logger.debug("Looking for SNPs near transcript %s (with an upstream promoter window size of %s)..." % (
        transcript_dict['transcript.id'], promoter_upstream_window_size))
    transcript_chrom = transcript_dict['assembly.ucsc_style_name'][3:]
    if transcript_chrom in snp_intervals:
        validation.validate_strand(transcript_dict['transcript.strand'])

        if transcript_dict['transcript.strand'] == '+':
            window_start = max(
                transcript_dict['transcript.start'] - 1 - promoter_upstream_window_size, 0)
            window_end = transcript_dict['transcript.end']
        else:
            window_start = transcript_dict['transcript.start']
            window_end = min(
                transcript_dict['transcript.end'], transcript_dict['assembly.sequence_length'])
        logger.debug('Looking for SNPS in region %s:%s-%s for transcript %s (strand %s)',
                     transcript_chrom, window_start, window_end, transcript_dict['transcript.id'],
                     transcript_dict['transcript.strand'])
        hits = set(snp_intervals[transcript_chrom][window_start:window_end])
    else:
        logger.debug('Transcript %s chromosome %s not found in SNP file.',
                     transcript_dict['transcript.id'], transcript_chrom)
        hits = set([])

    rs_ids = sorted(['rs' + interval.data for interval in hits])
    logger.debug('Found %s SNPs near transcript %s: %s',
                 len(hits), transcript_dict['transcript.id'], rs_ids)

    return rs_ids
