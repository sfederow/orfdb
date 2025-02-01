"""
This module provides tools for sequence analysis and manipulation. It includes functions for computing GC content, 
nucleotide probabilities, codon counts and frequencies, extracting sequence blocks, and more.
"""

import collections
import logging
import typing
from typing import Dict, List, Optional, Set, Union, Any, Tuple

import Bio
import Bio.Seq
import Bio.SeqIO

from orfdb.bigprot import misc, coordinates, validation, orf_classes, constants

logger = logging.getLogger(__name__)


UNIFORM_NUC_PROBS = {'A': 0.25,
                     'C': 0.25,
                     'G': 0.25,
                     'T': 0.25
                     }


def compute_gc_fraction(sequence: str) -> float:
    """
    Compute the GC content of a sequence.

    Parameters
    ----------
    sequence : str
        The sequence to compute the GC content for.

    Returns
    -------
    float
        The GC content of the sequence.
    """
    gc_count = 0
    non_standard_counts = collections.defaultdict(lambda: 0)

    for char in sequence:
        if char in {'G', 'C'}:
            gc_count += 1
        elif char not in {'A', 'T'}:
            non_standard_counts[char] += 1

    if non_standard_counts:
        logger.warning(
            'Warning -- found non-standard nucleotide characters: %s!', dict(non_standard_counts))

    gc_fraction = gc_count / \
        (len(sequence) - sum(non_standard_counts.values()))

    return gc_fraction


def compute_nuc_probs(sequence: str) -> dict[str, float]:
    """
    Compute the nucleotide probabilities of a sequence.

    Parameters
    ----------
    sequence : str
        The sequence to compute the nucleotide probabilities for.

    Returns
    -------
    dict[str, float]
        A dictionary mapping each nucleotide to its probability in the sequence.
    """
    gc_fraction = compute_gc_fraction(sequence)

    return generate_nuc_probs_from_gc(gc_fraction)


def generate_nuc_probs_from_gc(gc_fraction: float) -> dict[str, float]:
    """
    Generate nucleotide probabilities from a GC content.

    Parameters
    ----------
    gc_fraction : float
        The GC content to generate the nucleotide probabilities from.

    Returns
    -------
    dict[str, float]
        A dictionary mapping each nucleotide to its probability based on the GC content.
    """
    return {'A': (1 - gc_fraction) / 2,
            'C': gc_fraction / 2,
            'G': gc_fraction / 2,
            'T': (1 - gc_fraction) / 2
            }


def extract_seq_blocks(parent_seq: str, chrom_starts: list[int], block_sizes: list[int], source_coordinate_system: str = 'gff') -> list[str]:
    """
    Extract sequence blocks from a parent sequence given start positions and block sizes.
    Assumes gff-style intervals (1-based, fully-closed).

    Parameters
    ----------
    parent_seq : str
        The parent sequence to extract the blocks from.
    chrom_starts : list[int]
        The start positions of the blocks.
    block_sizes : list[int]
        The sizes of the blocks.
    source_coordinate_system : str
        The coordinate system of the start positions and block sizes.

    Returns
    -------
    list[str]
        A list of strings representing the extracted sequence blocks.
    """
    py_intervals = coordinates.convert_chrom_starts_and_block_sizes_to_intervals(chrom_starts=chrom_starts,
                                                                                 block_sizes=block_sizes,
                                                                                 source_type=source_coordinate_system,
                                                                                 destination_type='python')
    py_chrom_starts, py_block_ends = zip(*py_intervals)

    seq_blocks = []

    for chrom_start, block_end in zip(py_chrom_starts, py_block_ends):
        assert chrom_start >= 0, f'Start {chrom_start} must be greater than or equal to 0.'
        assert block_end <= len(
            parent_seq), f'End {block_end} must be less than {len(parent_seq)}'
        seq_blocks.append(parent_seq[chrom_start:block_end])

    return seq_blocks


def extract_sequence(parent_seq: str,
                     chrom_starts: List[int],
                     block_sizes: List[int],
                     strand: str,
                     source_coordinate_system: str = 'gff',
                     seq_record_id: str = '') -> Union[str, Bio.SeqIO.SeqRecord]:
    """
    Extract a sequence from a parent sequence given start positions, block sizes, and strand.

    Parameters
    ----------
    parent_seq : str
        The parent sequence to extract from.
    chrom_starts : list[int]
        The start positions of the blocks.
    block_sizes : list[int]
        The sizes of the blocks.
    source_coordinate_system : str
        The coordinate system of the start positions and block sizes.
    strand : str
        The strand ('+' or '-').
    seq_record_id : str, optional
        The sequence record ID.

    Returns
    -------
    Union[str, Bio.SeqIO.SeqRecord]
        The extracted sequence.
    """
    validation.validate_strand(strand)

    seq_blocks = extract_seq_blocks(
        parent_seq, chrom_starts, block_sizes, source_coordinate_system)
    try:
        full_seq = ''.join(seq_blocks)
    except TypeError:  # Assume the blocks are Bio.Seqs then
        full_seq = ''.join([str(block.seq) for block in seq_blocks])

    if strand == '-':
        full_seq = Bio.Seq.reverse_complement(full_seq)

    if seq_record_id:
        return Bio.SeqIO.SeqRecord(Bio.Seq.Seq(full_seq), id=seq_record_id)
    else:
        return full_seq


def extract_transcript_sequence_from_genome(transcript_dict: Dict, genome_dict: Dict[str, str], seq_record_id: str = '',
                                            source_coordinate_system: str = constants.TRANSCRIPT_COORDINATE_SYSTEM,
                                            accession_namespace: str = constants.DEFAULT_ACCESSION_NAMESPACE) -> Union[str, Bio.SeqIO.SeqRecord]:
    """
    Extract a transcript sequence from a genome given a transcript dictionary and a genome dictionary.

    Parameters
    ----------
    transcript_dict : dict
        The transcript dictionary.
    genome_dict : dict[str, str]
        The genome dictionary.
    seq_record_id : str, optional
        The sequence record ID.
    source_coordinate_system : str, optional
        The coordinate system of the start positions and block sizes (default is 'velia').
    accession_namespace : str, optional
        The accession namespace (default is 'genbank').

    Returns
    -------
    Union[str, Bio.SeqIO.SeqRecord]
        The extracted transcript sequence.
    """
    logger.debug('Extracting genomic sequence of transcript %s ...',
                 transcript_dict['transcript.id'])
    validation.validate_coordinate_type(source_coordinate_system)

    if transcript_dict['transcript.chrom_starts'] and transcript_dict['transcript.block_sizes']:
        chrom_starts, block_sizes = misc.extract_numeric_vector(
            transcript_dict['transcript.chrom_starts']), misc.extract_numeric_vector(
            transcript_dict['transcript.block_sizes'])
    else:  # if empty, make the whole transcript one block
        full_interval = coordinates.IntervalFactory.make_interval(source_coordinate_system)(start=int(transcript_dict['transcript.start']),
                                                                                            end=int(transcript_dict['transcript.end']))
        chrom_starts = [full_interval.python_start]
        block_sizes = [full_interval.python_size]
        source_coordinate_system = 'python'


    extracted_sequence = extract_sequence(
        parent_seq=genome_dict[transcript_dict[f'assembly.{accession_namespace}_accession']],
        chrom_starts=chrom_starts,
        block_sizes=block_sizes,
        strand=transcript_dict['transcript.strand'],
        source_coordinate_system=source_coordinate_system,
        seq_record_id=seq_record_id)

    return extracted_sequence


def extract_orf_sequence_from_genome(aliased_orf_dict: dict[str, typing.Any], genome_dict: dict[str, str], seq_record_id: str = '', source_coordinate_system: str = 'gff', accession_namespace: str = constants.DEFAULT_ACCESSION_NAMESPACE) -> str:
    """
    Extract an ORF sequence from a genome given an aliased ORF dictionary and a genome dictionary.
    Assumes gff-style intervals (1-based, fully-closed).

    Parameters
    ----------
    aliased_orf_dict : dict[str, typing.Any]
        The aliased ORF dictionary.
    genome_dict : dict[str, str]
        The genome dictionary.
    seq_record_id : str, optional
        The sequence record ID.
    source_coordinate_system : str, optional
        The coordinate system of the start positions and block sizes (default is 'gff').
    accession_namespace : str, optional
        The accession namespace (default is 'genbank').

    Returns
    -------
    str
        The extracted ORF sequence.
    """
    logger.debug('Extracting genomic sequence of ORF %s ...',
                 aliased_orf_dict['orf.orf_idx_str'])
    validation.validate_coordinate_type(source_coordinate_system)

    if aliased_orf_dict['orf.chrom_starts'] and aliased_orf_dict['orf.block_sizes']:
        chrom_starts, block_sizes = misc.extract_numeric_vector(
            aliased_orf_dict['orf.chrom_starts']), misc.extract_numeric_vector(
            aliased_orf_dict['orf.block_sizes'])
    else:  # if empty, make the whole ORF one block
        full_interval = coordinates.IntervalFactory.make_interval(source_coordinate_system)(start=int(aliased_orf_dict['transcript.start']),
                                                                                            end=int(aliased_orf_dict['transcript.end']))
        chrom_starts = [full_interval.python_start]
        block_sizes = [full_interval.python_size]
        source_coordinate_system = 'python'

    extracted_sequence = extract_sequence(parent_seq=genome_dict[aliased_orf_dict[f'assembly.{accession_namespace}_accession']],
                                          chrom_starts=chrom_starts,
                                          block_sizes=block_sizes,
                                          strand=aliased_orf_dict['orf.strand'],
                                          seq_record_id=seq_record_id,
                                          source_coordinate_system=source_coordinate_system)

    return extracted_sequence


def split_sequence_into_codons(nucleotide_sequence: str, frame: int = 0) -> list[str]:
    """
    Split a nucleotide sequence into codons given a frame.

    Parameters
    ----------
    nucleotide_sequence : str
        The nucleotide sequence to split.
    frame : int, optional
        The frame (default is 0).

    Returns
    -------
    list[str]
        A list of strings representing the codons.
    """
    validation.validate_frame(frame)

    codons = []
    for i in range((len(nucleotide_sequence) - frame) // 3):
        codon_start = i * 3 + frame
        codon_end = codon_start + 3
        codons.append(nucleotide_sequence[codon_start:codon_end])
    return codons


def extract_orf_start_context(orf_object: orf_classes.OrfBase,
                              parent_transcript_seq: str,
                              upstream_length: int,
                              downstream_length: int
                              ) -> tuple[str, str]:
    """
    Extract the start context of an ORF, which is the potential Kozak sequence starting 4 bp upstream and ending 1 bp downstream of the start codon.

    Parameters
    ----------
    orf_object : orf_classes.OrfBase
        The ORF object containing ORF information such as start position and nucleotide length.
    parent_transcript_seq : Bio.SeqIO.SeqRecord
        The parent transcript sequence from which the ORF start context will be extracted.

    Returns
    -------
    tuple[str, str]
        A tuple of strings representing the upstream and downstream sequences of the ORF start context.
    """
    assert orf_object.start_pos >= upstream_length, 'ORF %s start position must be greater or equal to %s for a valid Kozak sequence!' % (
        orf_object, upstream_length)
    assert orf_object.nuc_length >= downstream_length + \
        3, 'ORF %s nucleotide length must be greater or equal to %s for a valid Kozak sequence!' % (
            orf_object, downstream_length + 3)

    upstream_seq = parent_transcript_seq[orf_object.start_pos -
                                         upstream_length:orf_object.start_pos]
    downstream_seq = parent_transcript_seq[orf_object.start_pos +
                                           3:orf_object.start_pos+3 + downstream_length]
    assert len(upstream_seq) == upstream_length
    assert len(downstream_seq) == downstream_length

    return upstream_seq, downstream_seq


def find_internal_stop_codons(orf_dict: Dict, genome: Dict, stop_codons: Optional[List[str]] = None, accession_namespace: str = constants.DEFAULT_ACCESSION_NAMESPACE) -> List[
        Tuple[int, str]]:
    """Find internal stop codons in the ORF sequence.
    Args:
        orf_dict (Dict): Dictionary containing ORF information.
        genome (Dict): Dictionary containing genome information.
        stop_codons (Optional[List[str]]): List of stop codons. Defaults to ['TGA', 'TAA', 'TAG'].
    Returns:
        List[Tuple[int, str]]: List of tuples containing the index and codon of internal stop codons.
    """
    if stop_codons is None:
        stop_codons = constants.DEFAULT_STOP_CODONS
    stop_codons = set(stop_codons)
    orf_sequence = extract_orf_sequence_from_genome(aliased_orf_dict=orf_dict,
                                                    genome_dict=genome, source_coordinate_system='gff',
                                                    accession_namespace=accession_namespace)
    orf_phases = misc.extract_numeric_vector(orf_dict['orf.phases'])

    validation.validate_strand(orf_dict['orf.strand'])

    if orf_dict['orf.strand'] == '-':
        initial_orf_phase = orf_phases[-1]
    else:
        initial_orf_phase = orf_phases[0]

    validation.validate_phase(initial_orf_phase)

    if initial_orf_phase != 0:
        orf_sequence = orf_sequence[initial_orf_phase:]

    internal_stop_codons = []
    for codon_idx, codon in enumerate(split_sequence_into_codons(orf_sequence)[:-1]):
        if codon in stop_codons:
            internal_stop_codons.append((codon_idx, codon))

    return internal_stop_codons


def get_sequence_record(
    sequence: str,
    seq_record_id: str = ''
) -> Union[str, Bio.SeqIO.SeqRecord]:
    """Get sequence record from string."""
    # ... rest of the function ...


def get_sequence_from_fasta(
    fasta_fpath: str,
    accession_namespace: str = constants.DEFAULT_ACCESSION_NAMESPACE
) -> Union[str, Bio.SeqIO.SeqRecord]:
    """Get sequence from FASTA file."""
    # ... rest of the function ...
