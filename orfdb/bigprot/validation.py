"""Module for validating Open Reading Frames (ORFs).

This module provides functions for validating various aspects of ORFs such as coordinates, strand, frame, phase, etc.
"""
import logging
from typing import List, Dict, Optional, Set, Any

import numpy as np

from . import misc

VALID_COORDINATE_TYPES: Set[str] = {'python', 'gff', 'velia'}

logger = logging.getLogger(__name__)


class ExonValidationError(ValueError):
    """Custom exception for exon validation errors."""

    def __init__(self, message="There was an error with exon validation."):
        self.message = message
        super().__init__(self.message)


class ExonsNotContiguous(ExonValidationError):
    """Exception raised when exons are not contiguous."""

    def __init__(self, message="Exons are not contiguous."):
        self.message = message
        super().__init__(self.message)


class ExonsNotStartWithOne(ExonValidationError):
    """Exception raised when exon numbering does not start with 1."""

    def __init__(self, message="Exon numbering does not start with 1."):
        self.message = message
        super().__init__(self.message)


class ExonsTranscriptNotMatch(ExonValidationError):
    """Exception raised when an exon's transcript ID does not match the expected transcript ID."""

    def __init__(self, transcript_id=None, message="An exon's transcript ID does not match the expected transcript ID."):
        if transcript_id:
            message += f" Found transcript ID: {transcript_id}."
        self.message = message
        super().__init__(self.message)


class ExonCountNotRight(ExonValidationError):
    """Exception raised when the number of exons does not match the number of chromosome starts."""

    def __init__(self, expected_count=None, actual_count=None, message="The number of exons does not match the number of chromosome starts."):
        if expected_count and actual_count:
            message += f" Expected {expected_count}, but found {actual_count}."
        self.message = message
        super().__init__(self.message)


class TranscriptValidationError(ValueError):
    """Custom exception for transcript validation errors."""

    def __init__(self, message="There was an error with transcript validation."):
        self.message = message
        super().__init__(self.message)


class ChromStartsBlockSizesValidationError(ValueError):
    def __init__(self, message="There was an error with chrom starts and block sizes"):
        self.message = message
        super().__init__(self.message)


class ChromStartsBlockSizesOverlap(ChromStartsBlockSizesValidationError):
    """Exception raised when the blocks defined by chromosome starts and block sizes overlap."""

    def __init__(self, message="Chromosome starts and block sizes overlap."):
        self.message = message
        super().__init__(self.message)


class ChromStartsBlockSizesOutOfBounds(ChromStartsBlockSizesValidationError):
    """Exception raised when chromosome starts or block sizes are out of bounds."""

    def __init__(self, message="Chromosome starts or block sizes are out of bounds."):
        self.message = message
        super().__init__(self.message)


def validate_transcript(transcript_dict: Dict[str, Any]):
    transcript_start = transcript_dict['transcript.start'] - 1
    transcript_end = transcript_dict['transcript.end']

    try:
        validate_strand(transcript_dict['transcript.strand'])
        chrom_starts = misc.extract_numeric_vector(
            transcript_dict['transcript.chrom_starts'])
        block_sizes = misc.extract_numeric_vector(
            transcript_dict['transcript.block_sizes'])

        # Hacky but until we refactor the import structure we can't use coordinates module here
        chrom_starts = [cs - 1 for cs in chrom_starts]
        validate_chrom_starts_block_sizes(
            chrom_starts_python=chrom_starts, block_sizes_python=block_sizes)

    except (ChromStartsBlockSizesValidationError, AssertionError) as this_exception:
        raise TranscriptValidationError(this_exception.message)

    if chrom_starts[0] < transcript_start:
        raise TranscriptValidationError('First chrom start %s less than transcript start %s for transcript %s' % (
            chrom_starts[0], transcript_start, transcript_dict['transcript.id']))

    if chrom_starts[-1] + block_sizes[-1] > transcript_end:
        raise TranscriptValidationError('Last chrom start %s plus last block size %s gives %s which is greater than transcript end %s for transcript %s' % (
            chrom_starts[-1], block_sizes[-1], chrom_starts[-1] + block_sizes[-1], transcript_end, transcript_dict['transcript.id']))


def validate_exon_list(exon_list: List[Dict[str, Any]], parent_transcript: Dict[str, Any]) -> None:
    """Validates a list of exons for a given parent transcript.

    This function checks if all exons belong to the same transcript, are contiguous, start with exon number 1,
    and if the number of exons matches the number of chromosome starts in the parent transcript.

    Args:
        exon_list: A list of dictionaries, each representing an exon with keys like 'transcript_exon.transcript_id' and 'transcript_exon.exon_number'.
        parent_transcript: A dictionary representing the parent transcript with a key 'transcript.chrom_starts'.

    Raises:
        ExonsTranscriptNotMatch: If any exon's transcript ID does not match the expected transcript ID.
        ExonsNotContiguous: If exons are not contiguous.
        ExonsNotStartWithOne: If exon numbering does not start with 1.
        ExonCountNotRight: If the number of exons does not match the number of chromosome starts.
    """
    transcript_id = exon_list[0]['transcript_exon.transcript_id']
    logger.debug('Validating %s exons for transcript %s ...' %
                 (len(exon_list), transcript_id))

    for exon in exon_list:
        if exon['transcript_exon.transcript_id'] != transcript_id:
            raise ExonsTranscriptNotMatch(
                'Exon %s transcript ID does not match %s' % (exon, transcript_id))

    exon_nums = np.sort(
        np.array([exon['transcript_exon.exon_number'] for exon in exon_list]))

    if np.any(np.diff(exon_nums) != 1):
        raise ExonsNotContiguous(
            'Exons for transcript %s not contiguous!' % transcript_id)

    if exon_nums[0] != 1:
        raise ExonsNotStartWithOne(
            'Exons for transcript %s don\'t start with 1!' % transcript_id)

    if np.not_equal(np.diff(exon_nums), 1).sum() != 0:
        raise ExonsNotContiguous(
            'Exons for trranscript %s are not contiguous!' % transcript_id)

    num_chrom_starts = len(misc.extract_numeric_vector(
        parent_transcript['transcript.chrom_starts']))
    if num_chrom_starts != len(exon_list):
        raise ExonCountNotRight('Transcript %s has % exons but % chrom_starts!' % (
            transcript_id, len(exon_list), num_chrom_starts))


def validate_orf_coordinates(orf_dict: Dict, transcript_dict: Optional[Dict] = None):
    """Validate ORF coordinates.

    This function validates the ORF coordinates by comparing them with the inferred coordinates from the ORF's chrom_starts and block_sizes.
    If a transcript_dict is provided, it also validates the ORF coordinates against the parent transcript coordinates.

    Args:
        orf_dict (Dict): Dictionary containing ORF information.
        transcript_dict (Optional[Dict]): Dictionary containing parent transcript information. Defaults to None.

    Raises:
        AssertionError: If any of the ORF coordinates do not match the expected values.
    """
    if orf_dict['orf.chrom_starts'] and orf_dict['orf.block_sizes']:
        orf_chrom_starts = misc.extract_numeric_vector(orf_dict['orf.chrom_starts'],
                                                       remove_dangling_separators=True)
        orf_block_sizes = misc.extract_numeric_vector(orf_dict['orf.block_sizes'],
                                                      remove_dangling_separators=True)
        inferred_orf_start = orf_chrom_starts[0]
        inferred_orf_end = orf_chrom_starts[-1] + orf_block_sizes[-1]
        assert inferred_orf_start == orf_dict[
            'orf.start'], f'First ORF chrom_start {orf_chrom_starts[0]} not equal to annotated ORF start {orf_dict["orf.start"]}!'
        assert inferred_orf_end == orf_dict[
            'orf.end'], f'Last ORF chrom_start {orf_chrom_starts[-1]} plus last block size {orf_block_sizes[-1]} gives {inferred_orf_end} whch is not equal to annotated ORF start {orf_dict["orf.end"]}!'

    # If passed, validate against parent transcript coordinates
    if transcript_dict:
        assert orf_dict['orf.assembly_id'] == transcript_dict[
            'assembly.id'], f'ORF assembly ID {orf_dict["orf.assembly_id"]} does not match parent transcript assembly ID {transcript_dict["assembly.id"]}!'
        assert orf_dict['orf.strand'] == transcript_dict[
            'transcript.strand'], f'ORF strand {orf_dict["orf.strand"]} does not match parent transcriupt strand {transcript_dict["transcript.strand"]}!'
        assert orf_dict['orf.start'] >= transcript_dict[
            'transcript.start'], f'ORF start {orf_dict["orf.start"]} is upstream of parent transcript start {transcript_dict["transcript.start"]}!'
        assert orf_dict['orf.end'] <= transcript_dict[
            'transcript.end'], f'ORF end {orf_dict["orf.end"]} is upstream of parent transcript end {transcript_dict["transcript.end"]}!'


def validate_coordinate_type(coordinate_type: str):
    """Validate coordinate type.

    This function validates the coordinate type by checking if it is one of the valid coordinate types.

    Args:
        coordinate_type (str): The coordinate type to validate.

    Raises:
        ValueError: If the coordinate type is not valid.
    """
    if not coordinate_type in VALID_COORDINATE_TYPES:
        raise ValueError(
            f'Invalid coordinate type "{coordinate_type}. Valid types are: {",".join(VALID_COORDINATE_TYPES)}.')


def validate_strand(strand: str, valid_strands: Optional[Set[str]] = None):
    """Validate strand.

    This function validates the strand by checking if it is one of the valid strands.

    Args:
        strand (str): The strand to validate.
        valid_strands (Optional[Set[str]]): Set of valid strands. Defaults to {'+', '-'}.

    Raises:
        AssertionError: If the strand is not valid.
    """
    if valid_strands is None:
        valid_strands = {'+', '-'}
    assert strand in {
        '+', '-'}, f'Invalid strand value {strand}. Valid options are: {valid_strands}.'


def validate_frame(frame: int):
    """Validate frame.

    This function validates the frame by checking if it is one of the valid frames.

    Args:
        frame (int): The frame to validate.

    Raises:
        AssertionError: If the frame is not valid.
    """
    assert frame in {
        0, 1, 2}, f'Invalid frame value {frame}. Must be 0, 1 or 2.'


def validate_phase(phase: int):
    """Validate phase.

    This function validates the phase by checking if it is one of the valid phases.

    Args:
        phase (int): The phase to validate.

    Raises:
        AssertionError: If the phase is not valid.
    """
    assert phase in {
        0, 1, 2}, f'Invalid phase value {phase}. Must be 0, 1 or 2.'


def validate_phase_style(phase_style: str):
    assert phase_style in {
        'gencode', 'ensembl'}, f'Invalid phase style {phase_style}. Must be "gencode" or "ensembl".'


def validate_chrom_starts_block_sizes(chrom_starts_python: List[int], block_sizes_python: List[int]):
    """Validate chrom_starts and block_sizes.

    This function validates the chrom_starts and block_sizes by checking if they have the same number of elements.

    Args:
        chrom_starts (List[int]): Pythonic list of chrom_starts.
        block_sizes (List[int]): Pythonic list of block_sizes.

    Raises:
        AssertionError: If chrom_starts and block_sizes do not have the same number of elements.
    """
    try:
        assert len(chrom_starts_python) == len(
            block_sizes_python), f'Chrom starts and block sizes must have the same number of elements!'

        assert all(
            cs >= 0 for cs in chrom_starts_python), f'All chrom_starts must be non-negative!'
        assert all(
            bs >= 0 for bs in block_sizes_python), f'All block_sizes must be non-negative!'

    except AssertionError as ae:
        raise ChromStartsBlockSizesValidationError(ae.message)

    if not all([chrom_starts_python[i] >= chrom_starts_python[i-1] +
                block_sizes_python[i-1] for i in range(1, len(chrom_starts_python))]):
        raise ChromStartsBlockSizesOverlap(
            'Chrom starts must start after the end of the previous block!')
