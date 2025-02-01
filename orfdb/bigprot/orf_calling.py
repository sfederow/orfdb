"""
This module provides classes and functions for finding and annotating Open Reading Frames (ORFs) in a given sequence.
It includes classes for representing ORFs and sets of ORFs, and functions for finding ORFs in sequences and sequence sets.
"""

import logging
import multiprocessing
from datetime import datetime
from typing import List, Tuple, Iterable, Dict, Set, Collection, Any, Optional, Union

import Bio.SeqIO
import Bio.SeqRecord
import intervaltree

from orfdb.bigprot.orf_classes import OrfBase
from . import annotation
from . import bigwig_tracks
from . import constants
from . import misc
from . import motif_tools
from . import sequence_tools
from . import validation

logger = logging.getLogger(__name__)


def find_orfs_in_seq_strand(seq: Bio.SeqIO.SeqRecord,
                            strand: str,
                            min_codon_length: int = constants.DEFAULT_MIN_CODON_LENGTH,
                            max_codon_length: int = constants.DEFAULT_MAX_CODON_LENGTH,
                            start_codons: Collection[str] = constants.DEFAULT_START_CODONS,
                            stop_codons: Collection[str] = constants.DEFAULT_STOP_CODONS) -> Set[OrfBase]:
    """
    Find ORFs in a sequence on a given strand.

    :param seq: Sequence to search for ORFs.
    :param strand: Strand to search for ORFs.
    :param min_codon_length: Minimum length of an ORF in codons.
    :param max_codon_length: Maximum length of an ORF in codons.
    :param start_codons: List of start codons.
    :param stop_codons: List of stop codons.
    :return: Set of ORFs found in the sequence on the given strand.
    """
    start_codons = set(start_codons)
    stop_codons = set(stop_codons)

    validation.validate_strand(strand)
    assert min_codon_length >= 2, f'Got {min_codon_length} minimum length but ORFs must have a start and stop codon, so cannot have length < 2.'

    orf_class = OrfBase
    parent_sequence_id = seq.id

    this_seq = str(seq.seq)

    all_orfs = set([])

    for reading_frame in range(3):
        active_orfs = set([])

        for pos in range(reading_frame, len(this_seq) - 2, 3):
            this_codon = this_seq[pos:pos + 3]

            if this_codon in stop_codons:
                for this_orf in active_orfs:
                    this_orf.add_codon(this_codon)
                    if this_orf.codon_length >= min_codon_length:
                        if strand == '-':
                            # account for rev complementation
                            this_orf.start_pos = len(
                                this_seq) - this_orf.end_pos
                        all_orfs.add(this_orf)

                active_orfs = set([])
            else:
                # Stop extending ORFs as soon as they exceed the length limit (accounting for the fact that we need to end with a stop codon)
                discard_list = []
                for orf in active_orfs:
                    if orf.codon_length >= max_codon_length - 1:
                        discard_list.append(orf)
                    else:
                        orf.add_codon(this_codon)

                if len(discard_list):
                    active_orfs.difference_update(discard_list)

                # Start a new ORF if this is a start codon
                if this_codon in start_codons:
                    new_orf = orf_class(parent_sequence_id=parent_sequence_id,
                                        start_pos=pos,
                                        # we won't know the start position for minus strand ORFs until they're complete.
                                        strand=strand,
                                        codons=[this_codon])
                    active_orfs.add(new_orf)

    return all_orfs


def find_orfs_in_seq_set(seq_strand_iterable: Iterable[Tuple[Bio.SeqIO.SeqRecord, str]],
                         min_codon_length: int = constants.DEFAULT_MIN_CODON_LENGTH,
                         max_codon_length: int = constants.DEFAULT_MAX_CODON_LENGTH,
                         start_codons: List[str] = None,
                         stop_codons: List[str] = None
                         ) -> Dict[Union[str, None], Set[OrfBase]]:
    """
    Find ORFs in a set of sequences.

    :param seq_strand_iterable: Iterable of tuples, each containing a sequence and a strand.
    :param min_codon_length: Minimum length of an ORF in codons.
    :param max_codon_length: Maximum length of an ORF in codons.
    :param start_codons: List of start codons.
    :param stop_codons: List of stop codons.
    :return: Dictionary mapping sequence IDs to sets of ORFs found in the sequences.
    """
    if start_codons is None:
        start_codons = constants.DEFAULT_START_CODONS
    if stop_codons is None:
        stop_codons = constants.DEFAULT_STOP_CODONS

    orfs_by_seq_id = {}

    for this_seq, strand in seq_strand_iterable:
        if this_seq.id in orfs_by_seq_id:
            raise KeyError(
                f'SeqRecord with id {this_seq.id} has already been called!')

        orfs_by_seq_id[this_seq.id] = find_orfs_in_seq_strand(seq=this_seq,
                                                              strand=strand,
                                                              min_codon_length=min_codon_length,
                                                              max_codon_length=max_codon_length,
                                                              start_codons=start_codons,
                                                              stop_codons=stop_codons)

    return orfs_by_seq_id


def find_orfs_in_seq_set_mp(seq_strand_list: Collection[Tuple[Bio.SeqIO.SeqRecord, str]],
                            min_codon_length: int = constants.DEFAULT_MIN_CODON_LENGTH,
                            max_codon_length: int = constants.DEFAULT_MAX_CODON_LENGTH,
                            start_codons: Optional[List[str]] = None,
                            stop_codons: Optional[List[str]] = None,
                            num_processes: int = 0,
                            chunk_size: int = constants.DEFAULT_MP_CHUNK_SIZE) -> Dict[str, Set[OrfBase]]:
    """
    Finds Open Reading Frames (ORFs) in a collection of sequences using multiprocessing.

    :param seq_strand_list: A collection of tuples, each containing a sequence record and a strand.
    :param min_codon_length: The minimum length of an ORF in codons.
    :param max_codon_length: The maximum length of an ORF in codons.
    :param start_codons: A list of start codons. If None, defaults are used.
    :param stop_codons: A list of stop codons. If None, defaults are used.
    :param num_processes: The number of processes to use. If 0, uses the number of CPU cores.
    :param chunk_size: The size of chunks for dividing the work among processes.
    :return: A dictionary mapping sequence IDs to sets of ORFs found in the sequences.
    """
    if start_codons is None:
        start_codons = constants.DEFAULT_START_CODONS
    if stop_codons is None:
        stop_codons = constants.DEFAULT_STOP_CODONS
    if not num_processes:
        num_processes = multiprocessing.cpu_count()

    logger.info(
        'Finding ORFs in %d sequences using %d processes', len(seq_strand_list), num_processes)

    orfs_by_seq_id: Dict[int, Set[OrfBase]] = {}
    orf_pool = multiprocessing.Pool(num_processes)
    params = [(seq, strand, min_codon_length, max_codon_length, start_codons, stop_codons) for seq, strand in
              seq_strand_list]

    for seq_id, new_orfs in orf_pool.map(_find_orfs_mp_worker, params, chunksize=chunk_size):
        orfs_by_seq_id[seq_id] = new_orfs

    orf_pool.close()
    orf_pool.join()

    return orfs_by_seq_id


def _find_orfs_mp_worker(params: Tuple[Bio.SeqIO.SeqRecord, str, int, int, List[str], List[str]]) -> Tuple[
        str, Set[OrfBase]]:
    seq, strand, min_codon_length, max_codon_length, start_codons, stop_codons = params

    return (seq.id, find_orfs_in_seq_strand(seq=seq,
                                            strand=strand,
                                            min_codon_length=min_codon_length,
                                            max_codon_length=max_codon_length,
                                            start_codons=start_codons,
                                            stop_codons=stop_codons))


def call_and_annotate_orfs(transcripts_by_id: Dict[str, Dict[str, Any]],
                           genome: Dict[str, str],
                           output_prefix: str,
                           exons_by_transcript: Dict[str, Dict[int, List[Dict[str, Any]]]],
                           cdss_by_exon: Dict[str, List[Dict[str, Any]]],
                           skip_annotation: bool = False,
                           min_codon_length: int = constants.DEFAULT_MIN_CODON_LENGTH,
                           max_codon_length: int = constants.DEFAULT_MAX_CODON_LENGTH,
                           start_codons: List[str] = constants.DEFAULT_START_CODONS,
                           stop_codons: List[str] = constants.DEFAULT_STOP_CODONS,
                           phase_style: str = constants.DEFAULT_PHASE_STYLE,
                           accession_namespace: str = constants.DEFAULT_ACCESSION_NAMESPACE,
                           transcript_chunk_size: int = constants.DEFAULT_TRANSCRIPT_CHUNK_SIZE,
                           max_chunks: int = -1,
                           num_processes: int = constants.DEFAULT_NUM_PROCESSES) -> None:
    """
    This function identifies and annotates Open Reading Frames (ORFs) within a collection of transcripts, leveraging multiple processes for efficiency. It utilizes genomic sequences, PhyloCSF tracks for evolutionary conservation analysis, and SNP intervals for genetic variation insights. The function supports customization of ORF detection parameters, including codon length and start/stop codons, and allows for optional skipping of the annotation step.

    Parameters
    ----------
    transcripts_by_id : Dict[str, Dict[str, Any]]
        Maps transcript IDs to their associated metadata.
    genome : Dict[str, str]
        Associates chromosome names with their corresponding sequences.
    output_prefix : str
        Specifies the prefix for all generated output files.
    exons_by_transcript : Dict[str, Dict[int, List[Dict[str, Any]]]]
        Maps transcript IDs to their exons, facilitating detailed ORF annotation.
    cdss_by_exon : Dict[str, List[Dict[str, Any]]]
        Associates exon IDs with their coding sequences (CDSs), aiding in ORF annotation.
    skip_annotation : bool, optional
        If True, skips the annotation step to only identify ORFs. Defaults to False.
    min_codon_length : int
        Sets the minimum codon length for an ORF to be considered valid.
    max_codon_length : int
        Sets the maximum codon length for an ORF to be considered valid.
    start_codons : List[str]
        Defines the codons that can act as start codons for ORFs.
    stop_codons : List[str]
        Defines the codons that can act as stop codons, marking the end of ORFs.
    phase_style : str
        Determines the phase style to be used in ORF annotation.
    transcript_chunk_size : int
        Defines the size of transcript chunks for processing, optimizing performance.
    num_processes : int
        Specifies the number of processes to use for parallel computation.
    accession_namespace : str
        Defines the namespace for accession numbers in the genome sequences.

    Returns
    -------
    None
    """
    start_time = datetime.now()

    orf_table_output_fpath = output_prefix + '_orfs.csv.gz'
    supplemental_orf_table_output_fpath = output_prefix + '_orfs_supplemental.csv.gz'
    transcript_orf_table_output_fpath = output_prefix + '_transcript_orfs.csv.gz'
    cds_orf_table_output_fpath = output_prefix + '_cds_orfs.csv.gz'

    chunk_count = 0
    orf_count = 0
    logger.info(
        'Finding ORFs with codon lengths %d-%d with start codons %s and stop codons %s for %d transcripts in chunks of %d.', min_codon_length, max_codon_length, ",".join(start_codons), ",".join(stop_codons), len(transcripts_by_id), transcript_chunk_size)

  
    kozak_upstream_pssm, kozak_downstream_pssm = motif_tools.generate_upstream_downstream_kozak_87_pssm()

    transcript_keys = list(transcripts_by_id.keys())

    orf_table_writer = misc.LazyCsvWriter(
        output_fpath=orf_table_output_fpath, name='ORF table')
    supplemental_orf_table_writer = misc.LazyCsvWriter(
        output_fpath=supplemental_orf_table_output_fpath, name='Supplemental ORF table')
    transcript_orf_table_writer = misc.LazyCsvWriter(
        output_fpath=transcript_orf_table_output_fpath, name='Transcript ORF table')
    cds_orf_table_writer = misc.LazyCsvWriter(
        output_fpath=cds_orf_table_output_fpath, name='CDS ORF table')

    transcript_idx = 0
    current_orf_id = 1  # Orf IDs are 1-indexed in DB
    orf_ids_by_idx_str = {}
    cds_orf_linkages = set([])

    this_transcript_chunk = misc.get_dict_chunk(transcripts_by_id, transcript_keys[transcript_idx:min(
        transcript_idx + transcript_chunk_size, len(transcripts_by_id))])

    while this_transcript_chunk and (max_chunks < 0 or chunk_count < max_chunks):
        logger.info(
            'Processing chunk %d of %d with %d transcripts ...', chunk_count+1, len(transcripts_by_id) // transcript_chunk_size + 1, len(this_transcript_chunk))

        transcript_seq_strands_by_id = {
            str(transcript['transcript.id']): (sequence_tools.extract_transcript_sequence_from_genome(transcript,
                                                                                                      genome,
                                                                                                      source_coordinate_system=constants.TRANSCRIPT_COORDINATE_SYSTEM,
                                                                                                      seq_record_id=str(
                                                                                                          transcript[
                                                                                                              'transcript.id']),
                                                                                                      accession_namespace=accession_namespace
                                                                                                      ).upper(),
                                               transcript['transcript.strand'],
                                               ) for transcript in
            this_transcript_chunk.values()}

        this_chunk_orfs = find_orfs_in_seq_set_mp(transcript_seq_strands_by_id.values(),
                                                  min_codon_length=min_codon_length,
                                                  max_codon_length=max_codon_length,
                                                  start_codons=start_codons,
                                                  stop_codons=stop_codons,
                                                  num_processes=num_processes,
                                                  chunk_size=constants.DEFAULT_MP_CHUNK_SIZE)

        num_found_orfs = sum([len(this_transcript_orfs)
                              for this_transcript_orfs in this_chunk_orfs.values()])
        orf_count += num_found_orfs
        logger.info('Found %d (possibly redundant) ORFs ...', num_found_orfs)

        current_orf_id, this_orf_table_chunk, this_supplemental_orf_table_chunk, this_transcript_orf_table_chunk, this_cds_orf_table_chunk = annotation.annotate_orfs(
            current_orf_id=current_orf_id,
            transcript_orfs=this_chunk_orfs,
            transcripts_by_id=transcripts_by_id,
            transcript_seq_strands_by_id=transcript_seq_strands_by_id, genome=genome,
            orf_ids_by_idx_str=orf_ids_by_idx_str,
            exons_by_transcript=exons_by_transcript,
            cdss_by_exon=cdss_by_exon,
            cds_orf_linkages=cds_orf_linkages,
            kozak_upstream_pssm=kozak_upstream_pssm,
            kozak_downstream_pssm=kozak_downstream_pssm,
            phase_style=phase_style,
            accession_namespace=accession_namespace,
            skip_annotation=skip_annotation) 

        orf_table_writer.writerows(this_orf_table_chunk)
        supplemental_orf_table_writer.writerows(
            this_supplemental_orf_table_chunk)
        transcript_orf_table_writer.writerows(
            this_transcript_orf_table_chunk)
        cds_orf_table_writer.writerows(this_cds_orf_table_chunk)

        transcript_idx += transcript_chunk_size
        this_transcript_chunk = misc.get_dict_chunk(transcripts_by_id, transcript_keys[transcript_idx:min(
            transcript_idx + transcript_chunk_size, len(transcripts_by_id))])
        chunk_count += 1

    end_time = datetime.now()
    elapsed_time = end_time - start_time
    logger.info(
        'Done calling and annotating ORFs. Found %d ORFs, of which %d were unique, in %d transcripts.', orf_count, current_orf_id - 1, len(transcripts_by_id))
    logger.info('Job took %s.', elapsed_time)
