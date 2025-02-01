"""
This module contains tools for working with motifs in the context of Open Reading Frame (ORF) finding. 
It includes functionalities for converting count matrices to a format compatible with Biopython, 
and for computing motif scores. It also defines some constants related to the human genome and Kozak sequences.
"""

import logging
from collections import defaultdict
from typing import Tuple, Iterable, Dict, List

import Bio.Seq
import Bio.motifs
import Bio.motifs.matrix
import numpy as np
from tqdm import tqdm

from . import constants
from . import misc
from . import orf_classes
from . import sequence_tools

logger = logging.getLogger(__name__)

HG38_TRANSCRIPTOME_BACKGROUND: dict[str, float] = {
    'C': 0.23508065229347194,
    'T': 0.2565544499572134,
    'G': 0.23959576178514036,
    'A': 0.26876913596417434
}

HG38_GENOME_BACKGROUND: dict[str, float] = {
    'T': 0.29523216423206444,
    'A': 0.2943609111517895,
    'C': 0.2047778221112342,
    'G': 0.20562910250491193
}

# This was taken from BioStars (https://www.biostars.org/p/116533/) and thought to be confirmed in literature but cannot
# currently find literature source.
KOZAK_RAW_COUNT_ARRAY: np.ndarray = np.array([
    [25, 53, 15, 7],
    [61, 2, 36, 1],
    [27, 49, 13, 11],
    [15, 55, 21, 9],
    [97, 0, 0, 0],
    [0, 0, 0, 97],
    [0, 0, 97, 0],
    [23, 16, 46, 15]
])

# These are taken from https://doi.org/10.1093/nar/15.20.8125
KOZAK_87_UPSTREAM_MOTIF_PERCENTS = {'A': [23, 26, 25, 23, 19, 23, 17, 18, 25, 61, 27, 15],
                                    'C': [35, 35, 35, 26, 39, 37, 19, 39, 53, 2, 49, 55],
                                    'G': [26, 21, 22, 33, 23, 20, 44, 23, 15, 36, 13, 21],
                                    'T': [19, 18, 18, 18, 19, 20, 20, 20, 7, 1, 11, 9]}

KOZAK_87_DOWNSTREAM_MOTIF_PERCENTS = {'A': [23],
                                      'C': [16],
                                      'G': [46],
                                      'T': [15]}


def convert_count_matrix_to_biopython(count_matrix: np.ndarray[np.int64], alphabet=None) -> \
        dict[
            str, list[int]]:
    """
    Convert a count matrix to a format compatible with Biopython.

    Parameters
    ----------
    count_matrix : np.ndarray[np.int64]
        The count matrix to be converted. Each element in the matrix is an integer.
    alphabet : list[str], optional
        The alphabet of nucleotides, by default ['A', 'C', 'G', 'T']

    Returns
    -------
    dict[str, list[int]]
        The converted count matrix.
    """
    if alphabet is None:
        alphabet = ['A', 'C', 'G', 'T']
    assert count_matrix.shape[1] == len(alphabet)

    counts: dict[str, list[int]] = {nuc: [] for nuc in alphabet}
    for col_idx, nuc in enumerate(alphabet):
        for row_idx in range(count_matrix.shape[0]):
            counts[nuc].append(count_matrix[row_idx, col_idx])

    return counts


def convert_pcm_matrix_to_motif(pcm_matrix: np.ndarray[np.int64]) -> Bio.motifs.Motif:
    """
    Convert a position count matrix to a Biopython motif.

    Parameters
    ----------
    pcm_matrix : np.ndarray[int]
        The position count matrix to be converted.

    Returns
    -------
    Bio.motifs.Motif
        The converted Biopython motif.
    """
    return Bio.motifs.Motif(counts=convert_count_matrix_to_biopython(pcm_matrix))


def generate_pwm(motif: Bio.motifs.Motif,
                 pseudocount: float = constants.DEFAULT_MOTIF_PSEUDOCOUNT) -> Bio.motifs.matrix.PositionWeightMatrix:
    """
    Generate a position weight matrix (PWM) from a motif.

    Parameters
    ----------
    motif : Bio.motifs.Motif
        The motif to generate the PWM from.
    pseudocount : float, optional
        The pseudocount to be used, by default constants.DEFAULT_MOTIF_PSEUDOCOUNT

    Returns
    -------
    Bio.motifs.matrix.PositionWeightMatrix
        The generated PWM.
    """
    return motif.normalize(pseudocount)


def generate_pssm(pwm: Bio.motifs.matrix.PositionWeightMatrix,
                  background_model: dict[str, float]) -> Bio.motifs.matrix.PositionSpecificScoringMatrix:
    """
    Generate a position specific scoring matrix (PSSM) from a PWM and a background model.

    Parameters
    ----------
    pwm : Bio.motifs.matrix.PositionWeightMatrix
        The PWM to generate the PSSM from.
    background_model : dict[str, float]
        The background model to be used for generating the PSSM.

    Returns
    -------
    Bio.motifs.matrix.PositionSpecificScoringMatrix
        The generated PSSM.
    """
    return pwm.log_odds(background=background_model)


def count_nucs_over_sequences(sequences: Iterable[str], alphabet=None) -> \
        dict[str, int]:
    """
    Count the occurrences of nucleotides over a collection of sequences.

    Parameters
    ----------
    sequences : typing.Iterable[str]
        The collection of sequences to count nucleotides over.
    alphabet : typing.Iterable[str], optional
        The alphabet of nucleotides, by default ['A', 'C', 'G', 'T']

    Returns
    -------
    dict[str, int]
        The counts of nucleotides over the sequences.
    """
    if alphabet is None:
        alphabet = ['A', 'C', 'G', 'T']
    alphabet_set: set[str] = set(alphabet)
    all_nuc_counts: defaultdict = defaultdict(lambda: 0)
    for seq in tqdm(sequences):
        for nuc, count in misc.count(seq.upper()).items():
            if nuc in alphabet_set:
                all_nuc_counts[nuc] += count
    return all_nuc_counts


def convert_counts_to_freqs(count_dict: dict[str, int]) -> dict[str, float]:
    """
    Convert nucleotide counts to frequencies.

    Parameters
    ----------
    count_dict : dict[str, int]
        The dictionary containing nucleotide counts.

    Returns
    -------
    dict[str, float]
        The dictionary containing nucleotide frequencies.
    """
    total_count: int = sum(count_dict.values())
    return {key: count / total_count for key, count in count_dict.items()}


def compute_motif_score(sequence: str, pssm: Bio.motifs.matrix.PositionSpecificScoringMatrix) -> float:
    """
    Compute the motif score for a given sequence using a position specific scoring matrix (PSSM).

    Parameters
    ----------
    sequence : str
        The sequence for which the motif score is to be computed.
    pssm : Bio.motifs.matrix.PositionSpecificScoringMatrix
        The position specific scoring matrix (PSSM) to be used for computing the motif score.

    Returns
    -------
    float
        The computed motif score.
    """
    assert len(sequence) == pssm.length, 'Length of passed sequence %s doesn\'t match length of PSSM %s!' % (
        len(sequence), pssm.length)

    try:
        pos, score = next(pssm.search(
            sequence, threshold=float('-inf'), both=False))
    except StopIteration:
        logger.warning(
            'No match for PSSM %s in sequence %s!', pssm, sequence)
        return np.NaN
    else:
        assert pos == 0
        return score


def compute_kozak_score(
        orf_object: orf_classes.OrfBase,
        parent_transcript_seq: str,
        kozak_upstream_pssm: Bio.motifs.matrix.PositionSpecificScoringMatrix,
        kozak_downstream_pssm: Bio.motifs.matrix.PositionSpecificScoringMatrix
) -> float:
    """
    Compute the Kozak score for a given ORF using upstream and downstream Kozak sequences.
    That is, we omit the "ATG" portion of the kozak sequence, since we want to allow other start codons.
    So we take the sum of the scores for the upstream and downstream portions of the Kozak sequence that flank
    the start codon.

    Parameters
    ----------
    orf_object : orf_classes.OrfBase
        The ORF object containing ORF information.
    parent_transcript_seq : str
        The sequence of the parent transcript from which the Kozak sequence will be extracted.
    kozak_upstream_pssm : Bio.motifs.matrix.PositionSpecificScoringMatrix
        The position specific scoring matrix (PSSM) for the upstream Kozak sequence.
    kozak_downstream_pssm : Bio.motifs.matrix.PositionSpecificScoringMatrix
        The position specific scoring matrix (PSSM) for the downstream Kozak sequence.

    Returns
    -------
    float
        The Kozak sequence motif score for the given ORF.
    """
    logger.debug('Computing Kozak sequence score for %s ...',
                 orf_object)
    if orf_object.start_pos < kozak_upstream_pssm.length or orf_object.nuc_length < kozak_downstream_pssm.length + 3:
        logger.debug(
            'ORF %s does not have sufficient flanking sequence to compute a Kozak motif score.', orf_object)
        return np.nan

    upstream_kozak_seq, downstream_kozak_seq = sequence_tools.extract_orf_start_context(
        orf_object=orf_object, parent_transcript_seq=parent_transcript_seq, upstream_length=kozak_upstream_pssm.length,
         downstream_length=kozak_downstream_pssm.length)

    upstream_score: float = compute_motif_score(
        sequence=upstream_kozak_seq, pssm=kozak_upstream_pssm)
    logger.debug('Computed upstream Kozak score of %s', upstream_score)

    downstream_score: float = compute_motif_score(
        sequence=downstream_kozak_seq, pssm=kozak_downstream_pssm)
    logger.debug('Computed downstream Kozak score of %s', downstream_score)

    overall_score: float = upstream_score + downstream_score
    logger.debug('Overall Kozak score: %s', overall_score)

    return overall_score


def generate_upstream_downstream_kozak_84_pssm(kozak_pcm: np.ndarray[np.int64] = KOZAK_RAW_COUNT_ARRAY,
                                               background_model: Dict[str, float] = HG38_TRANSCRIPTOME_BACKGROUND,
                                               pseudocount: float = constants.DEFAULT_MOTIF_PSEUDOCOUNT) -> Tuple[
        Bio.motifs.matrix.PositionSpecificScoringMatrix, Bio.motifs.matrix.PositionSpecificScoringMatrix]:
    """
    Generate upstream and downstream Kozak PSSMs from a position count matrix (PCM).

    Parameters
    ----------
    kozak_pcm : np.ndarray, optional
        The position count matrix for the Kozak sequence, by default KOZAK_RAW_COUNT_ARRAY
    background_model : dict[str, float], optional
        The background model for nucleotide frequencies, by default HG38_TRANSCRIPTOME_BACKGROUND
    pseudocount : float, optional
        The pseudocount to be added to the counts before normalization, by default constants.DEFAULT_MOTIF_PSEUDOCOUNT

    Returns
    -------
    Tuple[Bio.motifs.matrix.PositionSpecificScoringMatrix, Bio.motifs.matrix.PositionSpecificScoringMatrix]
        A tuple containing the upstream and downstream Kozak PSSMs.
    """
    if background_model is None:
        background_model = HG38_TRANSCRIPTOME_BACKGROUND
    kozak_motif_upstream = Bio.motifs.Motif(
        counts=convert_count_matrix_to_biopython(kozak_pcm[:4, :]))
    kozak_motif_downstream = Bio.motifs.Motif(
        counts=convert_count_matrix_to_biopython(kozak_pcm[7:8, :]))

    kozak_pssm_upstream = kozak_motif_upstream.counts.normalize(
        pseudocount).log_odds(background=background_model)
    kozak_pssm_downstream = kozak_motif_downstream.counts.normalize(
        pseudocount).log_odds(background=background_model)

    return kozak_pssm_upstream, kozak_pssm_downstream


def generate_upstream_downstream_kozak_87_pssm(kozak_upstream_pcm: Dict[str, List[int]] = KOZAK_87_UPSTREAM_MOTIF_PERCENTS,
                                               kozak_downstream_pcm: Dict[str, List[int]
                                                                          ] = KOZAK_87_UPSTREAM_MOTIF_PERCENTS,
                                               background_model=None,
                                               pseudocount: float = constants.DEFAULT_MOTIF_PSEUDOCOUNT) -> Tuple[
        Bio.motifs.matrix.PositionSpecificScoringMatrix, Bio.motifs.matrix.PositionSpecificScoringMatrix]:
    """
    Generate upstream and downstream Kozak PSSMs from a position count matrix (PCM).

    Parameters
    ----------
    kozak_pcm : np.ndarray, optional
        The position count matrix for the Kozak sequence, by default KOZAK_RAW_COUNT_ARRAY
    background_model : dict[str, float], optional
        The background model for nucleotide frequencies, by default HG38_TRANSCRIPTOME_BACKGROUND
    pseudocount : float, optional
        The pseudocount to be added to the counts before normalization, by default constants.DEFAULT_MOTIF_PSEUDOCOUNT

    Returns
    -------
    Tuple[Bio.motifs.matrix.PositionSpecificScoringMatrix, Bio.motifs.matrix.PositionSpecificScoringMatrix]
        A tuple containing the upstream and downstream Kozak PSSMs.
    """
    if background_model is None:
        background_model = HG38_TRANSCRIPTOME_BACKGROUND
    kozak_motif_upstream = Bio.motifs.Motif(
        counts=kozak_upstream_pcm)
    kozak_motif_downstream = Bio.motifs.Motif(
        counts=kozak_downstream_pcm)

    kozak_pssm_upstream = kozak_motif_upstream.counts.normalize(
        pseudocount).log_odds(background=background_model)
    kozak_pssm_downstream = kozak_motif_downstream.counts.normalize(
        pseudocount).log_odds(background=background_model)

    return kozak_pssm_upstream, kozak_pssm_downstream
