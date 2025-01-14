"""
This module contains functions for computing null codon probabilities, 
number of exclusive offsets, expected starts, and for generating and shuffling sequences.
"""

import itertools
import logging

import numpy as np

from . import sequence_tools

logger = logging.getLogger(__name__)


def compute_null_codon_probability(codon: str, nuc_probs: dict[str, float]) -> float:
    """
    Compute the probability of a codon given nucleotide probabilities.

    Parameters
    ----------
    codon : str
        The codon to compute the probability for.
    nuc_probs : dict[str, float]
        The probabilities of the nucleotides.

    Returns
    -------
    float
        The probability of the codon.
    """
    prob = 1
    for nuc in codon:
        prob *= nuc_probs[nuc]
    return prob


def compute_num_exclusive_offsets(codon_a: str, codon_b: str, start_offset: int = -2, end_offset: int = 2) -> int:
    """
    Compute the number of exclusive offsets between two codons.

    Parameters
    ----------
    codon_a : str
        The first codon.
    codon_b : str
        The second codon.
    start_offset : int, optional
        The start offset, by default -2.
    end_offset : int, optional
        The end offset, by default 2.

    Returns
    -------
    int
        The number of exclusive offsets.
    """
    exclusive_offset_count = 0

    for offset in range(start_offset, end_offset + 1):
        if offset == 0:
            sub_a = codon_a
            sub_b = codon_b
        elif offset > 0:
            sub_a = codon_a[offset:]
            sub_b = codon_b[:-offset]
        else:
            sub_a = codon_a[:offset]
            sub_b = codon_b[-offset:]

        if sub_a != sub_b:
            exclusive_offset_count += 1
        # print(offset, sub_a, sub_b, sub_a == sub_b)
    return exclusive_offset_count


def compute_expected_starts(sequence_length: int, nuc_probs: dict[str, float], start_codons: list[str]) -> float:
    """
    Compute the expected number of start codons in a sequence.

    Parameters
    ----------
    sequence_length : int
        The length of the sequence.
    nuc_probs : dict[str, float]
        The probabilities of the nucleotides.
    start_codons : list[str]
        The list of start codons.

    Returns
    -------
    float
        The expected number of start codons.
    """
    all_codon_probabilities = {start_codon: compute_null_codon_probability(start_codon, nuc_probs=nuc_probs) for
                               start_codon in start_codons}

    expected_starts = 0.0

    for start_codon in start_codons:
        # print(start_codon, compute_null_codon_likelihood(start_codon))
        expected_starts += all_codon_probabilities[start_codon] * \
                           sequence_length

    # print(expected_starts)
    for codon_1, codon_2 in itertools.combinations(start_codons, 2):
        if sequence_length < 7:
            for pos in range(0, sequence_length - 2):
                min_left_offset = max(-pos, -2)
                max_right_offset = min(sequence_length - pos - 3, 2)
                expected_starts -= all_codon_probabilities[codon_1] * all_codon_probabilities[
                    codon_2] * compute_num_exclusive_offsets(codon_1, codon_2, start_offset=min_left_offset,
                                                             end_offset=max_right_offset) / (
                                           max_right_offset - min_left_offset + 1)
        else:
            # Do the left edge
            max_right_offset = 2
            for min_left_offset in range(-2, 1):
                expected_starts -= all_codon_probabilities[codon_1] * all_codon_probabilities[
                    codon_2] * compute_num_exclusive_offsets(codon_1, codon_2, start_offset=min_left_offset,
                                                             end_offset=max_right_offset) / (
                                           max_right_offset - min_left_offset + 1)
            # Do the right edge
            min_left_offset = -2
            for max_right_offset in range(0, 3):
                expected_starts -= all_codon_probabilities[codon_1] * all_codon_probabilities[
                    codon_2] * compute_num_exclusive_offsets(codon_1, codon_2, start_offset=min_left_offset,
                                                             end_offset=max_right_offset) / (
                                           max_right_offset - min_left_offset + 1)
            # Account for the middle
            expected_starts -= all_codon_probabilities[codon_1] * all_codon_probabilities[
                codon_2] * compute_num_exclusive_offsets(codon_1, codon_2, start_offset=-2,
                                                         end_offset=2) * (sequence_length - 7) / (
                                       max_right_offset - min_left_offset + 1)

    # print(expected_starts)
    return expected_starts


def compute_expected_starts_from_sequence(sequence: str,
                                          start_codons: list[str]) -> float:
    """
    Compute the expected number of start codons from a sequence.

    Parameters
    ----------
    sequence : str
        The sequence to compute the expected number of start codons from.
    start_codons : list[str]
        The list of start codons.

    Returns
    -------
    float
        The expected number of start codons.
    """
    this_nuc_probs = sequence_tools.compute_nuc_probs(sequence=sequence)
    # print(this_nuc_probs)
    return compute_expected_starts(sequence_length=len(sequence),
                                   nuc_probs=this_nuc_probs,
                                   start_codons=start_codons)


def generate_random_sequence(seq_length: int, nuc_probs: dict[str, float], seed: int = 0) -> str:
    """
    Generate a random sequence of a given length.

    Parameters
    ----------
    seq_length : int
        The length of the sequence to generate.
    nuc_probs : dict[str, float]
        The probabilities of the nucleotides.
    seed : int, optional
        The seed for the random number generator, by default 0.

    Returns
    -------
    str
        The generated sequence.
    """
    if seed:
        np.random.seed(seed)

    return ''.join(np.random.choice(a=list(nuc_probs.keys()), size=seq_length, p=list(nuc_probs.values())))


def shuffle_string(sequence: str) -> str:
    """
    Shuffle a sequence.

    Parameters
    ----------
    sequence : str
        The sequence to shuffle.

    Returns
    -------
    str
        The shuffled sequence.
    """
    shuff_seq = np.array(list(sequence))
    np.random.shuffle(shuff_seq)
    return ''.join(shuff_seq)


def shuffle_sequence(sequence: str, num_shuffles: int = 1, seed: int = 0) -> list[str]:
    """
    Shuffle a sequence a number of times.

    Parameters
    ----------
    sequence : str
        The sequence to shuffle.
    num_shuffles : int, optional
        The number of times to shuffle the sequence, by default 1.
    seed : int, optional
        The seed for the random number generator, by default 0.

    Returns
    -------
    list[str]
        The list of shuffled sequences.
    """
    if seed:
        np.random.seed(seed)

    shuffle_seq = list(sequence)
    shuffled_seqs = []

    for _ in range(num_shuffles):
        np.random.shuffle(shuffle_seq)
        shuffled_seqs.append(''.join(shuffle_seq))

    return shuffled_seqs


def shuffled_sequence_generator(sequence: str, seed: int = 0) -> str:
    """
    Generate shuffled sequences from a sequence.

    Parameters
    ----------
    sequence : str
        The sequence to generate shuffled sequences from.
    seed : int, optional
        The seed for the random number generator, by default 0.

    Yields
    ------
    str
        The shuffled sequences.
    """
    if seed:
        np.random.seed(seed)

    shuffle_seq = list(sequence)
    while True:
        np.random.shuffle(shuffle_seq)
        yield ''.join(shuffle_seq)
