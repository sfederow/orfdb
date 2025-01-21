"""
This module handles operations related to BigWig tracks. It includes a class BwTracks that encapsulates the functionality 
for handling BigWig files, including loading the tracks, determining shared chromosomes, and other related operations.
"""
import glob
import itertools
import logging
from tqdm import tqdm
from typing import Dict, List, Set, Tuple

import numpy as np
import pyBigWig

import orfdb.bigprot.misc
from . import coordinates
from . import validation

logger = logging.getLogger(__name__)


class BwTracks:
    """
    A class to handle operations related to BigWig tracks, modified to use shared memory for multiprocessing.

    Attributes:
        trackset_name (str): The name of the trackset.
        bw_file_glob (str): The glob pattern to match BigWig files.
        _strand_position (int): The position of the strand in the file path. Default is -5.
        _frame_position (int): The position of the frame in the file path. Default is -4.
        bw_fpaths (list): List of file paths matched by the glob pattern.
        bw_tracks (list): List of BigWig tracks.
        bw_tracks_by_strand_frame (dict): Dictionary of BigWig tracks indexed by strand and frame.
        bw_vectors_by_strand_frame (dict): Dictionary of BigWig vectors indexed by strand and frame.
        shared_chroms (set): Set of chromosomes shared by all BigWig tracks.
        chrom_lengths (dict): Dictionary of chromosome lengths.
    """

    def __init__(self, bw_file_glob: str, trackset_name: str, load_to_ram: bool = False, strand_position: int = -5,
                 frame_position: int = -4, make_16bit=False):
        """
        Initializes the BwTracks object.

        Args:
            bw_file_glob (str): The glob pattern to match BigWig files.
            trackset_name (str): The name of the trackset.
            strand_position (int, optional): The position of the strand in the file path. Default is -5.
            frame_position (int, optional): The position of the frame in the file path. Default is -4.
        """
        self.trackset_name: str = trackset_name
        self.bw_file_glob: str = bw_file_glob
        self.load_to_ram: bool = load_to_ram
        self._strand_position: int = strand_position
        self._frame_position: int = frame_position
        self._make_16bit: bool = make_16bit

        self.bw_fpaths: List[str] = glob.glob(bw_file_glob)

        self.bw_tracks: List[pyBigWig.BigWigFile] = []

        self.bw_tracks_by_strand_frame: Dict[Tuple[str, int],
                                             pyBigWig.BigWigFile] = {}
        self.bw_vectors_by_strand_frame: Dict[Tuple[str, int],
                                              Dict[str, np.ndarray[np.float32]]] = {}
        self.shared_chroms: Set[str] = set([])
        self.chrom_lengths: Dict[str, int] = {}

        self._load_bw_tracks()

    def _check_completness(self):
        STRANDS = ('+', '-')
        FRAMES = (0, 1, 2)
        for strand, frame in itertools.product(STRANDS, FRAMES):
            if not (strand, frame) in self.bw_tracks_by_strand_frame:
                raise ValueError(
                    "Missing track for strand %s and frame %s in dataset %s" % (strand, frame, self.trackset_name))

    def _load_bw_tracks(self) -> None:
        """
        Loads BigWig tracks from the file paths.
        """
        logger.debug('Loading BigWig tracks for phylocsf tracks %s from %s ...',
                     self.trackset_name, self.bw_file_glob)

        if not self.bw_fpaths:
            raise ValueError(
                'No BigWig files found matching the glob pattern %s.' % self.bw_file_glob)

        for fpath in sorted(self.bw_fpaths):
            strand = fpath[self._strand_position]
            frame = int(fpath[self._frame_position]) - 1
            logger.debug(
                f'Loading bigwig track data {self.trackset_name} for strand {strand}, frame {frame - 1} from {fpath} ...')

            this_track = pyBigWig.open(fpath)
            if strand not in self.bw_tracks_by_strand_frame:
                self.bw_tracks_by_strand_frame[(strand, frame)] = {}
            self.bw_tracks_by_strand_frame[(strand, frame)] = this_track
            self.bw_tracks.append(this_track)

        self._determine_shared_chroms_bw_tracks()
        self._check_completness()

        if self.load_to_ram:
            for strand_frame in self.bw_tracks_by_strand_frame:
                self.bw_vectors_by_strand_frame[strand_frame] = {}
                chrom_length_dict = self.bw_tracks_by_strand_frame[(
                    strand, frame)].chroms()
                logger.info('Loading contents of bigwig vector for strand %s, frame %s into memory ...' % (
                    strand_frame[0], strand_frame[1]))
                for chrom in self.shared_chroms:
                    length = chrom_length_dict[chrom]
                    self.bw_vectors_by_strand_frame[strand_frame][chrom] = self.bw_tracks_by_strand_frame[
                        strand_frame].values(chrom,
                                             0, length, numpy=True)
                    if self._make_16bit:
                        self.bw_vectors_by_strand_frame[strand_frame][chrom] = \
                            self.bw_vectors_by_strand_frame[strand_frame][chrom].astype(
                                np.float16)

    def _determine_shared_chroms_bw_tracks(self) -> None:
        """
        Determines the chromosomes shared by all BigWig tracks.
        """
        logger.debug(
            'Determining set of shared chromosomes between BigWig tracks ...')
        first_pass = True
        for bw_track in self.bw_tracks:
            if first_pass:
                self.shared_chroms.update(
                    {chrom: None for chrom in bw_track.chroms()})
                first_pass = False
            else:
                self.shared_chroms.intersection_update(bw_track.chroms())
        all_chrom_lengths = bw_track.chroms()
        self.chrom_lengths = {
            chrom: all_chrom_lengths[chrom] for chrom in self.shared_chroms}
        logger.debug('Found %d shared chromosomes.', len(self.shared_chroms))

    def extract_bw_vector(self, strand: str, frame: int, chrom: str, py_start: int, py_end: int) -> np.ndarray[
            np.float32]:
        """
        Extracts the BigWig values for a given strand, frame, chromosome, and start and end positions.

        Note: Phylocsf enumerates their frames as {1,2,3} and we internally use {0,1,2}.

        Args:
            strand (str): The strand.
            frame (int): The frame.
            chrom (str): The chromosome.
            py_start (int): The start position.
            py_end (int): The end position.

        Returns:
            list: The BigWig values.
        """
        if self.load_to_ram:
            bw_vector = self.bw_vectors_by_strand_frame[(
                strand, frame)][chrom][py_start:py_end]
        else:
            bw_vector = self.bw_tracks_by_strand_frame[(
                strand, frame)].values(chrom, py_start, py_end, numpy=True)

        logger.debug('Extracted bigwig vector of length %s from dataset %s, frame %s, chrom %s, range %s-%s.',
                     len(bw_vector), self.trackset_name, frame, chrom, py_start, py_end)

        return bw_vector

    def extract_orf_bw_vector(self, orf_dict: Dict[str, str], parent_transcript_dict: Dict[str, str]) -> np.ndarray:
        """
        Extracts the BigWig values for a given ORF and parent transcript.

        Args:
            orf_dict (dict): The ORF dictionary.
            parent_transcript_dict (dict): The parent transcript dictionary.

        Returns:
            np.ndarray: The BigWig values.
        """
        logger.debug('Extracting bigwig vector for ORF %s (%s:%d-%d) from dataset %s ...',
                     orf_dict['orf.orf_idx_str'], parent_transcript_dict['assembly.ucsc_style_name'],
                      orf_dict['orf.start'], orf_dict['orf.end'], self.trackset_name)

        strand = orf_dict['orf.strand']
        validation.validate_strand(strand)

        chrom_length = int(parent_transcript_dict['assembly.sequence_length'])
        chrom = parent_transcript_dict['assembly.ucsc_style_name']

        gff_chrom_starts = orf_finding.misc.extract_numeric_vector(orf_dict['orf.chrom_starts'],
                                                                   remove_dangling_separators=True)
        gff_block_sizes = orf_finding.misc.extract_numeric_vector(orf_dict['orf.block_sizes'],
                                                                  remove_dangling_separators=True)
        validation.validate_chrom_starts_block_sizes(
            gff_chrom_starts, gff_block_sizes)

        orf_cds_intervals = coordinates.compute_orf_cds_intervals(orf_dict)

        py_chrom_starts, py_block_sizes = coordinates.convert_chrom_starts_and_block_sizes(gff_chrom_starts,
                                                                                           gff_block_sizes,
                                                                                           'gff',
                                                                                           'python')

        _, computed_frames = coordinates.compute_cds_phases_and_frames(py_chrom_starts, py_block_sizes, strand,
                                                                       chrom_length)
        logger.debug('Computed frames %s for ORF %s',
                     computed_frames, orf_dict['orf.orf_idx_str'])
        bw_value_chunks = []
        assert len(computed_frames) == len(orf_cds_intervals)
        
        for frame, (py_cds_start, py_cds_end) in zip(computed_frames, orf_cds_intervals):
            if chrom in self.shared_chroms:
                bw_value_chunks.append(self.extract_bw_vector(strand=strand,
                                                              frame=frame,
                                                              chrom=chrom,
                                                              py_start=py_cds_start,
                                                              py_end=py_cds_end))

            else:
                logger.debug('Chrom %s not found in dataset %s.',
                             chrom, self.trackset_name)
                bw_value_chunks.append(
                    [np.nan] * (py_cds_end - py_cds_start))

        return np.concatenate(bw_value_chunks)

    def compute_summary_stats(self, value_array: np.ndarray) -> Dict[str, float]:
        """
        Computes summary statistics for a given array of BigWig values.

        Args:
            value_array (np.ndarray): The array of BigWig values.

        Returns:
            dict: The summary statistics.
        """
        assert len(value_array) > 0, 'Value array has 0 length!'

        return {self.trackset_name + '_mean': value_array.mean(),
                self.trackset_name + '_median': np.median(value_array),
                self.trackset_name + '_std': value_array.std(),
                self.trackset_name + '_max': value_array.max(),
                self.trackset_name + '_min': value_array.min(),
                self.trackset_name + '_coverage': 1 - np.isnan(value_array).sum() / len(value_array)
                }


def bigwigger(bw_tracks_a: 'BwTracks', bw_tracks_b: 'BwTracks', operator_method_name: str) -> Dict[Tuple[str, int], Dict[str, np.ndarray]]:
    """
    Processes two sets of BigWig tracks to apply a specified operation on shared chromosomes and strand/frames.

    Args:
        bw_tracks_a (BwTracks): The first set of BigWig tracks.
        bw_tracks_b (BwTracks): The second set of BigWig tracks.
        operator_method_name (str): The name of the operation method to apply (e.g., 'add', 'subtract').

    Returns:
        Dict[Tuple[str, int], Dict[str, np.ndarray]]: A dictionary with keys as strand/frame tuples and values as dictionaries
        mapping chromosome names to the result of applying the operation on the chromosome vectors.
    """
    shared_strand_frames = sorted(set(bw_tracks_a.bw_tracks_by_strand_frame.keys(
    )).intersection(bw_tracks_b.bw_tracks_by_strand_frame.keys()))
    if not shared_strand_frames:
        logger.warning('No shared strand/frames!')
        return {}

    shared_chroms = sorted(
        bw_tracks_a.shared_chroms.intersection(bw_tracks_b.shared_chroms))
    if not shared_chroms:
        logger.warning('No shared chromosomes!')
        return {}

    result_tracks = {}
    for strand_frame in shared_strand_frames:
        logger.info(
            f'Processing strand {strand_frame[0]} frame {strand_frame[1]}')
        result_tracks[strand_frame] = {}
        chrom_length_dict = bw_tracks_a.bw_tracks_by_strand_frame[strand_frame].chroms(
        )
        for chrom in tqdm(shared_chroms):
            chrom_vector_a = bw_tracks_a.bw_tracks_by_strand_frame[strand_frame].values(
                chrom, 0, chrom_length_dict[chrom], numpy=True)
            chrom_vector_b = bw_tracks_b.bw_tracks_by_strand_frame[strand_frame].values(
                chrom, 0, chrom_length_dict[chrom], numpy=True)

            result_tracks[strand_frame][chrom] = getattr(
                chrom_vector_a, operator_method_name)(chrom_vector_b)
    return result_tracks


def write_non_nan_to_bigwig_vectorized(np_array, bw_file, chrom, start=0, step=1):
    """
    Write non-NaN portions of a NumPy array to an open bigWig file in fixedStep format
    using vectorized operations for improved performance.

    Parameters:
    - np_array: NumPy array containing the values to write.
    - bw_file: An open bigWig file object from pyBigWig.
    - chrom: The chromosome name as a string.
    - start: The start position (0-based) in the genome for the first element in np_array.
    - step: The step size between elements in np_array, assuming uniform spacing.

    Example usage:

    bw = pyBigWig.open("example.bw", "w")
    bw.addHeader([("chr1", 249250621)])  # Adjust this header based on your actual data

    Example NumPy array with some NaN values:
    np_array = np.array([0.1, np.nan, 0.2, 0.3, np.nan, np.nan, 0.4])

    Write the non-NaN data to the bigWig file, assuming uniform spacing:
    write_non_nan_to_bigwig_vectorized(np_array, bw, "chr1", start=100000, step=1000)

    bw.close()
    """
    # Mask to identify non-NaN elements
    valid_mask = ~np.isnan(np_array)

    # Find the start indices of contiguous non-NaN regions
    edges = np.diff(valid_mask.astype(int))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0] + 1

    # Handle case where the array starts with a non-NaN value
    if valid_mask[0]:
        starts = np.insert(starts, 0, 0)
    # Handle case where the array ends with a non-NaN value
    if valid_mask[-1]:
        ends = np.append(ends, valid_mask.size)

    # Calculate the genomic start positions and lengths of non-NaN regions
    genomic_starts = start + starts * step
    lengths = (ends - starts) * step

    # Write each non-NaN region to the bigWig file in fixedStep format
    for genomic_start, length in zip(genomic_starts, lengths):
        # Calculate the end position of the current region
        genomic_end = genomic_start + length
        # Extract the values for the current non-NaN region
        values = np_array[starts[0]:ends[0]][valid_mask[starts[0]:ends[0]]]
        # Add the entry to the bigWig file
        bw_file.addEntries([chrom] * len(values), [genomic_start] * len(values),
                           ends=[genomic_end] * len(values), values=values.tolist(), validate=False)
