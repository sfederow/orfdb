"""
This module provides classes and functions for handling intervals in different coordinate systems.
It includes an abstract base class for intervals, concrete classes for intervals in Python, GFF, and Velia coordinates,
a factory class for creating intervals of different types, and functions for converting between different coordinate systems.
"""
import logging
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Type, Union, Any

from . import misc, orf_classes, validation, constants

logger = logging.getLogger(__name__)


class BaseInterval(ABC):
    """
    Abstract base class for interval types.
    """

    @abstractmethod
    def __init__(self, start: int, size: int, end: int) -> None:
        """
        Initialize an interval.

        Args:
            start (int): Start of the interval.
            size (int): Size of the interval.
            end (int): End of the interval.
        """
        self._start = start
        self._size = size
        self._end = end

    def __repr__(self) -> str:
        return f'Interval [{self._start},{self._end}), size {self._size}'

    def _check_internal_consistency(self) -> None:
        """
        Check if the interval is internally consistent.
        """
        assert self._start >= 0, 'Start value %s cannot be negative!' % self._start
        assert self._end >= self._start, 'End value %s must be greater than start value %s!' % (
            self._start, self._end)
        assert self._end - \
            self._start == self._size, f'Interval {self} is not internally consistent!'

    @property
    def python_start(self) -> int:
        """
        Get the start of the interval in Python coordinates.

        Returns:
            int: Start of the interval.
        """
        return self._start

    @property
    def python_end(self) -> int:
        """
        Get the end of the interval in Python coordinates.

        Returns:
            int: End of the interval.
        """
        return self._end

    @property
    def python_size(self) -> int:
        """
        Get the size of the interval in Python coordinates.

        Returns:
            int: Size of the interval.
        """
        return self._size

    @property
    def gff_start(self) -> int:
        """
        Get the start of the interval in GFF coordinates.

        Returns:
            int: Start of the interval.
        """
        return self._start + 1

    @property
    def gff_size(self) -> int:
        """
        Get the size of the interval in GFF coordinates.

        Returns:
            int: Size of the interval.
        """
        return self._size

    @property
    def gff_end(self) -> int:
        """
        Get the end of the interval in GFF coordinates.

        Returns:
            int: End of the interval.
        """
        return self._end

    @property
    def velia_start(self) -> int:
        """
        Get the start of the interval in Velia coordinates.

        Returns:
            int: Start of the interval.
        """
        return self._start + 1

    @property
    def velia_end(self) -> int:
        """
        Get the end of the interval in Velia coordinates.

        Returns:
            int: End of the interval.
        """
        return self._end

    @property
    def velia_size(self) -> int:
        """
        Get the size of the interval in Velia coordinates.

        Returns:
            int: Size of the interval.
        """
        return self._size - 1

    @property
    def python_interval(self) -> Tuple[int, int]:
        """
        Get the interval in Python coordinates.

        Returns:
            tuple: Start and end of the interval.
        """
        return self.python_start, self.python_end

    @property
    def python_block(self) -> Tuple[int, int]:
        """
        Get the block in Python coordinates.

        Returns:
            tuple: Start and size of the block.
        """
        return self.python_start, self.python_size

    @property
    def gff_interval(self) -> Tuple[int, int]:
        """
        Get the interval in GFF coordinates.

        Returns:
            tuple: Start and end of the interval.
        """
        return self.gff_start, self.gff_end

    @property
    def gff_block(self) -> Tuple[int, int]:
        """
        Get the block in GFF coordinates.

        Returns:
            tuple: Start and size of the block.
        """
        return self.gff_start, self.gff_size

    @property
    def velia_interval(self) -> Tuple[int, int]:
        """
        Get the interval in Velia coordinates.

        Returns:
            tuple: Start and end of the interval.
        """
        return self.velia_start, self.velia_end

    @property
    def velia_block(self) -> Tuple[int, int]:
        """
        Get the block in Velia coordinates.

        Returns:
            tuple: Start and size of the block.
        """
        return self.velia_start, self.velia_size

    def convert_to_specified_block_type(self, block_type: str) -> Tuple[int, int]:
        """
        Convert the interval to a specified block type.

        Args:
            block_type (str): The block type to convert to.

        Returns:
            tuple: Start and size of the block in the specified coordinates.
        """
        validation.validate_coordinate_type(block_type)

        if block_type == 'gff':
            return self.gff_block
        elif block_type == 'velia':
            return self.velia_block
        else:
            return self.python_block

    def convert_to_specified_interval_type(self, interval_type: str) -> Tuple[int, int]:
        """
        Convert the interval to a specified interval type.

        Args:
            interval_type (str): The interval type to convert to.

        Returns:
            tuple: Start and end of the interval in the specified coordinates.
        """
        validation.validate_coordinate_type(interval_type)

        if interval_type == 'gff':
            return self.gff_interval
        elif interval_type == 'velia':
            return self.velia_interval
        else:
            return self.python_interval


class PythonInterval(BaseInterval):
    """
    Class for intervals in Python coordinates.
    """

    def __init__(self, start: int, end: int = -1, size: int = -1) -> None:
        """
        Initialize a Python interval.

        Args:
            start (int): Start of the interval.
            end (int, optional): End of the interval. Defaults to -1.
            size (int, optional): Size of the interval. Defaults to -1.

        Raises:
            ValueError: If neither end nor size is provided.
        """
        # Only 1 argument passed
        if end == -1 and size == -1:
            raise ValueError(
                'Must provide either an end coordinate or a size!')

        self._start = start

        if end != -1:
            self._end = end
        else:
            self._end = start + size

        if size != -1:
            self._size = size
        else:
            self._size = end - start

        self._check_internal_consistency()


class GffInterval(BaseInterval):
    """
    Class for intervals in GFF coordinates.
    """

    def __init__(self, start: int, end: int = -1, size: int = -1) -> None:
        """
        Initialize a GFF interval.

        Args:
            start (int): Start of the interval.
            end (int, optional): End of the interval. Defaults to -1.
            size (int, optional): Size of the interval. Defaults to -1.

        Raises:
            ValueError: If neither end nor size is provided.
        """
        # Only 1 argument passed
        if end == -1 and size == -1:
            raise ValueError(
                'Must provide either an end coordinate or a size!')

        self._start = start - 1  # Convert to 1-based indexing

        if end != -1:
            self._end = end
        else:
            self._end = start + size - 1  # Account for fully closed intervals

        if size != -1:
            self._size = size
        else:
            self._size = end - start + 1  # Account for fully closed intervals

        self._check_internal_consistency()


class VeliaInterval(BaseInterval):
    """
    Class for intervals in Velia coordinates.
    """

    def __init__(self, start: int, end: int = -1, size: int = -1) -> None:
        """
        Initialize a Velia interval.

        Args:
            start (int): Start of the interval.
            end (int, optional): End of the interval. Defaults to -1.
            size (int, optional): Size of the interval. Defaults to -1.

        Raises:
            ValueError: If neither end nor size is provided.
        """
        # Only 1 argument passed
        if end == -1 and size == -1:
            raise ValueError(
                'Must provide either an end coordinate or a size!')

        self._start = start - 1  # Convert to 1-based indexing

        # Seems like an odd way of doing this, but it allows us to always use the two existing (and internally consistent)
        # attributes to compute the third, which simplifies reasoning and debugging.
        if end != -1:
            self._end = end
        if size != -1:
            self._size = size + 1

        if end == -1:
            self._end = self._start + self._size
        if size == -1:
            self._size = self._end - self._start

        self._check_internal_consistency()


class IntervalFactory:
    """
    Factory class for creating intervals of different types.
    """

    @staticmethod
    def make_interval(interval_type: str) -> Type[Union[VeliaInterval, GffInterval, PythonInterval]]:
        """
        Create an interval of the specified type.

        Args:
            interval_type (str): The type of interval to create.

        Returns:
            Type[Union[VeliaInterval, GffInterval, PythonInterval]]: An interval of the specified type.
        """
        validation.validate_coordinate_type(interval_type)
        if interval_type == 'gff':
            interval_class = GffInterval
        elif interval_type == 'velia':
            interval_class = VeliaInterval
        else:
            interval_class = PythonInterval

        return interval_class


def convert_chrom_starts_and_block_sizes_to_interval_objects(chrom_starts: List[int], block_sizes: List[int],
                                                             interval_type: str) -> List[BaseInterval]:
    """
    Convert chromosome starts and block sizes to interval objects.

    Args:
        chrom_starts (List[int]): List of chromosome start positions.
        block_sizes (List[int]): List of block sizes.
        interval_type (str): The type of interval to create.

    Returns:
        List[BaseInterval]: List of interval objects.
    """
    validation.validate_coordinate_type(interval_type)
    assert len(chrom_starts) == len(block_sizes)

    interval_class = IntervalFactory.make_interval(interval_type)

    return [interval_class(start=cs, size=bs) for cs, bs in zip(chrom_starts, block_sizes)]


def convert_chrom_starts_and_block_sizes(chrom_starts: List[int], block_sizes: List[int], source_type: str,
                                         destination_type: str) -> Tuple[List[int], List[int]]:
    """
    Convert chromosome starts and block sizes from one coordinate system to another.

    Args:
        chrom_starts (List[int]): List of chromosome start positions.
        block_sizes (List[int]): List of block sizes.
        source_type (str): The type of the source coordinate system.
        destination_type (str): The type of the destination coordinate system.

    Returns:
        Tuple[List[int], List[int]]: Converted chromosome starts and block sizes.
    """
    validation.validate_coordinate_type(source_type)
    validation.validate_coordinate_type(destination_type)
    assert len(chrom_starts) == len(block_sizes)

    intervals = convert_chrom_starts_and_block_sizes_to_interval_objects(chrom_starts=chrom_starts,
                                                                         block_sizes=block_sizes,
                                                                         interval_type=source_type)
    converted_chrom_starts: List[int] = []
    converted_block_sizes: List[int] = []
    for this_interval in intervals:
        cs, bs = this_interval.convert_to_specified_block_type(
            destination_type)
        converted_chrom_starts.append(cs)
        converted_block_sizes.append(bs)

    return converted_chrom_starts, converted_block_sizes


def convert_chrom_starts_and_block_sizes_to_intervals(chrom_starts: List[int], block_sizes: List[int], source_type: str,
                                                      destination_type: str) -> List[Tuple[int, int]]:
    """
    Convert chromosome starts and block sizes to intervals.

    Args:
        chrom_starts (List[int]): List of chromosome start positions.
        block_sizes (List[int]): List of block sizes.
        source_type (str): The type of the source coordinate system.
        destination_type (str): The type of the destination coordinate system.

    Returns:
        List[Tuple[int, int]]: List of intervals.
    """
    validation.validate_coordinate_type(source_type)
    validation.validate_coordinate_type(destination_type)

    intervals = convert_chrom_starts_and_block_sizes_to_interval_objects(chrom_starts=chrom_starts,
                                                                         block_sizes=block_sizes,
                                                                         interval_type=source_type)

    return [this_interval.convert_to_specified_interval_type(destination_type) for this_interval in intervals]


def compute_orf_genomic_start(orf_object: orf_classes.OrfBase,
                              parent_transcript_dict: Dict[str, Union[str, int]]) -> int:
    """
    Computes the genomic start of an ORF object.

    This function assumes 'orf_object' has python coordinates and `parent_transcript_dict` has Velia coordinates.

    Args:
        orf_object: The ORF object.
        parent_transcript_dict: The parent transcript dictionary.

    Returns:
        The genomic start of the ORF object.
    """
    validation.validate_strand(orf_object.strand)
    gff_transcript_chrom_starts = misc.extract_numeric_vector(
        parent_transcript_dict['transcript.chrom_starts'])
    gff_transcript_block_sizes = misc.extract_numeric_vector(
        parent_transcript_dict['transcript.block_sizes'])

    py_transcript_chrom_starts, py_transcript_block_sizes = convert_chrom_starts_and_block_sizes(
        chrom_starts=gff_transcript_chrom_starts,
        block_sizes=gff_transcript_block_sizes,
        source_type=constants.TRANSCRIPT_COORDINATE_SYSTEM,
        destination_type='python')

    transcript_length = sum(py_transcript_block_sizes)
    accumulated_transcript_length = 0
    transcript_block_idx = 0

    while accumulated_transcript_length < transcript_length:
        transcript_chrom_start, transcript_block_size = py_transcript_chrom_starts[transcript_block_idx], \
            py_transcript_block_sizes[transcript_block_idx]
        transcript_offset = orf_object.start_pos - accumulated_transcript_length
        accumulated_transcript_length += transcript_block_size

        if orf_object.start_pos < accumulated_transcript_length:
            # orf must start within this transcript block
            orf_genomic_start = transcript_chrom_start + transcript_offset
            return orf_genomic_start

        transcript_block_idx += 1

    raise ValueError(f"We've run out of transcript length before we've exhausted the ORF length. Something's not "
                     f"right with the coordinates in transcript {parent_transcript_dict} or ORF {orf_object}.")


def compute_orf_cds_intervals(orf_dict: Dict[str, Union[str, int]]) -> List[Tuple[int, int]]:
    """
    Computes the CDS intervals of an ORF.

    Args:
        orf_dict: The ORF dictionary.

    Returns:
        The exon intervals of the ORF.
    """
    return convert_chrom_starts_and_block_sizes_to_intervals(
        misc.extract_numeric_vector(
            orf_dict['orf.chrom_starts'], remove_dangling_separators=True),
        misc.extract_numeric_vector(
            orf_dict['orf.block_sizes'], remove_dangling_separators=True),
        source_type='gff', destination_type='python')


def compute_orf_relative_start(gff_orf_genomic_start: int, parent_transcript_dict: Dict[str, Union[str, int]]) -> int:
    """
    Computes the start position of an ORF relative to the transcript sequence, given its genomic start position.

    This function assumes `orf_genomic_start` is in python coordinates and `parent_transcript_dict` has Velia coordinates.

    Args:
        gff_orf_genomic_start: The genomic start position of the ORF.
        parent_transcript_dict: The parent transcript dictionary.

    Returns:
        The start position of the ORF relative to the transcript sequence.
    """
    gff_transcript_chrom_starts = misc.extract_numeric_vector(
        parent_transcript_dict['transcript.chrom_starts'])
    gff_transcript_block_sizes = misc.extract_numeric_vector(
        parent_transcript_dict['transcript.block_sizes'])

    py_transcript_chrom_starts, py_transcript_block_sizes = convert_chrom_starts_and_block_sizes(
        chrom_starts=gff_transcript_chrom_starts,
        block_sizes=gff_transcript_block_sizes,
        source_type='gff',
        destination_type='python')

    py_genomic_start, _ = GffInterval(
        gff_orf_genomic_start, gff_orf_genomic_start + 1).python_interval

    accumulated_transcript_length = 0
    orf_relative_start = 0

    for chrom_start, block_size in zip(py_transcript_chrom_starts, py_transcript_block_sizes):
        if chrom_start <= py_genomic_start < chrom_start + block_size:
            # ORF starts within this block
            orf_relative_start += py_genomic_start - \
                chrom_start + accumulated_transcript_length
            return orf_relative_start
        accumulated_transcript_length += block_size
        orf_relative_start += accumulated_transcript_length

    raise ValueError(
        f"Could not compute the ORF's relative start position. Check the coordinates in transcript {parent_transcript_dict} or ORF genomic start {gff_orf_genomic_start}.")


def compute_orf_chrom_starts_block_sizes_exons(orf_object: orf_classes.OrfBase, orf_genomic_start: int,
                                               parent_transcript_dict: Dict[str, Union[str, int]]) -> Tuple[
        List[int], List[int], List[int]]:
    """
    Computes the chromosomal start positions, block sizes, and exon indices of an ORF.

    This function assumes `orf_object` has python coordinates and `parent_transcript_dict` has velia coordinates.
    It calculates the chromosomal start positions and block sizes for the ORF based on its genomic start position
    and the parent transcript's exon structure. Additionally, it identifies the indices of exons that overlap with
    the ORF.

    Args:
        orf_object: The ORF object, which should contain python coordinates.
        orf_genomic_start: The genomic start position of the ORF in python coordinates.
        parent_transcript_dict: The parent transcript dictionary containing velia coordinates.

    Returns:
        A tuple containing three lists:
        - The first list contains the chromosomal start positions of the ORF.
        - The second list contains the block sizes of the ORF.
        - The third list contains the indices of the exons that overlap with the ORF.
    """
    validation.validate_strand(orf_object.strand)
    orf_length = orf_object.end_pos - orf_object.start_pos

    gff_transcript_chrom_starts = misc.extract_numeric_vector(
        parent_transcript_dict['transcript.chrom_starts'])
    gff_transcript_block_sizes = misc.extract_numeric_vector(
        parent_transcript_dict['transcript.block_sizes'])

    py_transcript_chrom_starts, py_exon_sizes = convert_chrom_starts_and_block_sizes(
        chrom_starts=gff_transcript_chrom_starts,
        block_sizes=gff_transcript_block_sizes,
        source_type='gff',
        destination_type='python')

    computed_orf_chrom_starts, computed_orf_block_sizes, exons = [], [], []

    remaining_orf_length = orf_length
    current_exon_idx = 0

    while remaining_orf_length:
        exon_start, exon_size = py_transcript_chrom_starts[current_exon_idx], \
            py_exon_sizes[current_exon_idx]
        exon_end = exon_start + exon_size

        if exon_end > orf_genomic_start:
            # This exon overlaps the ORF
            this_orf_chrom_start = max(
                orf_genomic_start, exon_start)
            this_orf_block_size = min(exon_end - max(this_orf_chrom_start, exon_start),
                                      remaining_orf_length)

            computed_orf_chrom_starts.append(this_orf_chrom_start)
            computed_orf_block_sizes.append(this_orf_block_size)
            exons.append(current_exon_idx)
            remaining_orf_length -= this_orf_block_size

        assert current_exon_idx < len(
            py_exon_sizes), 'We have run out of exons before running out of ORF length!'
        current_exon_idx += 1

    assert sum(computed_orf_block_sizes) == orf_length

    return computed_orf_chrom_starts, computed_orf_block_sizes, exons


def compute_cds_phases_and_frames(py_chrom_starts: List[int], py_block_sizes: List[int], strand: str,
                                  chrom_length: int, phase_style: str = constants.DEFAULT_PHASE_STYLE) -> Tuple[List[int], List[int]]:
    """
    Computes the end phase of each block (assuming the blocks are CDSs) in the same manner as Ensembl.

    Note that the Ensembl GFF does not contain the end phase of the last CDS (because it's nonsensical) but
    we currently output it anyway.

    Args:
        py_chrom_starts: The chromosomal start positions.
        py_block_sizes: The block sizes.
        strand: The strand.
        chrom_length: The length of the chromosome.

    Returns:
        A tuple containing the phases and reading frames.
    """
    validation.validate_chrom_starts_block_sizes(chrom_starts_python=py_chrom_starts,
                                                 block_sizes_python=py_block_sizes)
    validation.validate_phase_style(phase_style)

    block_coords = convert_chrom_starts_and_block_sizes_to_intervals(chrom_starts=py_chrom_starts,
                                                                     block_sizes=py_block_sizes,
                                                                     source_type='python',
                                                                     destination_type='python')

    reading_frames = []
    phases = []
    # Start with the first CDS having a phase 0

    validation.validate_strand(strand)
    if strand == '+':
        intron_lengths = [block_coords[i + 1][0] - block_coords[i][1]
                          for i in range(len(block_coords) - 1)]

        for i, block in enumerate(block_coords):
            if i == 0:
                phases.append(0)
                reading_frames.append(block[0] % 3)
            else:
                intron_length = intron_lengths[i - 1]
                reading_frames.append(
                    (reading_frames[i - 1] + intron_length) % 3)
                phases.append((py_block_sizes[i-1] + phases[i-1]) % 3)

    elif strand == '-':
        block_coords = block_coords[::-1]
        py_block_sizes = py_block_sizes[::-1]

        intron_lengths = [block_coords[i][0] - block_coords[i + 1][1]
                          for i in range(len(block_coords) - 1)]

        for i, block in enumerate(block_coords):
            if i == 0:
                phases.append(0)
                reading_frames.append((chrom_length - block[1]) % 3)
            else:
                intron_length = intron_lengths[i - 1]
                reading_frames.append(
                    (reading_frames[i - 1] + intron_length) % 3)
                phases.append((py_block_sizes[i-1] + phases[i-1]) % 3)
        phases.reverse()
        reading_frames.reverse()

    if phase_style == 'gencode':
        phases = [(3-phase) % 3 for phase in phases]

    return phases, reading_frames
