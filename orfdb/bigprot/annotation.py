"""
This module contains functions for computing and validating Open Reading Frames (ORFs).

It includes functions to compute genomic start of an ORF, exon intervals, chromosomal start positions and block sizes,
end phase of each block, ORF value, ORF index, and PhyloCSF statistics for an ORF. It also includes functions to
extract and encode numeric vectors, and to populate the VELIA database fields for an ORF.
"""
import hashlib
import logging
from typing import List, Tuple, Union, Dict, Set, Any

import Bio.motifs.matrix
import intervaltree

from . import bigwig_tracks
from . import constants
from . import coordinates
from . import gwas_tools
from . import misc
from . import motif_tools
from . import orf_classes
from . import sequence_tools

logger = logging.getLogger(__name__)


def compute_orf_val(assembly_id: str,
                    start: int,
                    end: int,
                    strand: str,
                    block_sizes: str,
                    chrom_starts: str,
                    ) -> str:
    """
    Computes the ORF index string.

    Args:
        assembly_id: The assembly ID.
        chrom_starts: The chromosomal start positions.
        start: The overall start position.
        end: The overall end position.
        strand: The strand.
        block_sizes: The block sizes.

    Returns:
        A string uniquely representing this ORF..
    """
    return f'{assembly_id}_{start}_{end}_{strand}_{chrom_starts}_{block_sizes}'


def compute_orf_idx(orf_val: str) -> str:
    """
    Computes the ORF index.

    Args:
        orf_val: The ORF value.

    Returns:
        The ORF index.
    """
    return hashlib.sha3_512(orf_val.encode('utf-8')).hexdigest()


def populate_velia_db_fields_orf(orf_object: orf_classes.OrfBase,
                                 parent_transcript_dict: Dict[str, Union[str, int]],
                                 orf_id: str = '',
                                 velia_id_value: int = -1,
                                 phase_style=constants.DEFAULT_PHASE_STYLE,
                                 include_table_name: bool = True) -> Tuple[
        Dict[str, Union[int, str]], List[int], List[int]]:
    """
    Populates the VELIA database fields for an ORF.

    Args:
        orf_object: The ORF object.
        parent_transcript_dict: The parent transcript dictionary.
        orf_id: The ORF ID. Defaults to an empty string.
        velia_id_value: The VELIA ID value. Defaults to -1.
        include_table_name: Whether to include the table name. Defaults to True.

    Returns:
        A tuple containing:
            orf_attributes: the populated VELIA database fields for the ORF as a dictionary, 
            exon_indices: A list of integers indicating which of the enumerated exons associated with this parent transcript are overlapped by this ORF
            orf_cds_phases: A list of integers indicating the phase of each CDS in this ORF, one per exon that was overlapped by the ORF. 
            orf_cds_frames: A list of integers indicating the reading frame of each CDS in this ORF, one per exon that was overlapped by the ORF
    """
    logger.debug(
        'Populating Velia database fields for ORF %s ...', str(orf_object))
    assembly_id = parent_transcript_dict['assembly.id']

    py_genomic_start = coordinates.compute_orf_genomic_start(orf_object=orf_object,
                                                             parent_transcript_dict=parent_transcript_dict)
    logger.debug('Computed pythonic genomic start %s.', py_genomic_start)

    py_chrom_starts, py_block_sizes, exon_indices = coordinates.compute_orf_chrom_starts_block_sizes_exons(
        orf_object=orf_object,
        orf_genomic_start=py_genomic_start,
        parent_transcript_dict=parent_transcript_dict)
    logger.debug('Computed pythonic chrom starts %s and block sizes %s.',
                 py_chrom_starts, py_block_sizes)

    gff_chrom_starts, gff_block_sizes = coordinates.convert_chrom_starts_and_block_sizes(
        chrom_starts=py_chrom_starts,
        block_sizes=py_block_sizes,
        source_type='python',
        destination_type='gff')
    logger.debug('Computed GFF-style chrom starts %s and block sizes %s.',
                 gff_chrom_starts, gff_block_sizes)

    py_genomic_end = py_chrom_starts[-1] + py_block_sizes[-1]
    logger.debug('Computed pythonic genomic end %s', py_genomic_end)

    gff_genomic_start, gff_genomic_end = coordinates.PythonInterval(start=py_genomic_start,
                                                                    end=py_genomic_end).gff_interval
    logger.debug('Computed GFF-style start-end: %s-%s',
                 gff_genomic_start, gff_genomic_end)

    assert gff_genomic_start == gff_chrom_starts[0], 'Incorrect ORF start for ORF %s: %s' % (
        orf_object, gff_genomic_start)

    orf_cds_phases, orf_cds_frames = coordinates.compute_cds_phases_and_frames(py_chrom_starts=py_chrom_starts,
                                                                    py_block_sizes=py_block_sizes,
                                                                    strand=orf_object.strand,
                                                                    chrom_length=parent_transcript_dict[
                                                                        'assembly.sequence_length'],
                                                                    phase_style=phase_style)
    logger.debug('Got CDS phases %s and frames %s for ORF %s',
                 orf_cds_phases, orf_cds_frames, orf_object)

    gff_chrom_starts_string = misc.encode_numeric_vector(gff_chrom_starts)
    gff_block_sizes_string = misc.encode_numeric_vector(gff_block_sizes)

    orf_idx_string = compute_orf_val(start=gff_genomic_start,
                                     end=gff_genomic_end,
                                     strand=orf_object.strand,
                                     assembly_id=assembly_id,
                                     block_sizes=gff_block_sizes_string,
                                     chrom_starts=gff_chrom_starts_string,
                                     )

    orf_attributes = {'id': orf_id,
                      'start': gff_genomic_start,
                      'end': gff_genomic_end,
                      'strand': orf_object.strand,
                      'assembly_id': assembly_id,
                      'block_sizes': gff_block_sizes_string,
                      'chrom_starts': gff_chrom_starts_string,
                      'phases': misc.encode_numeric_vector(orf_cds_phases),
                      'exon_frames': misc.encode_numeric_vector(orf_cds_frames),
                      'orf_idx': compute_orf_idx(orf_idx_string),
                      'orf_idx_str': orf_idx_string,
                      'secondary_orf_id': '',
                      'benchling_id': '',
                      'aa_seq': orf_object.aa_sequence,
                      'nt_seq': orf_object.nuc_sequence,
                      'ensembl_protein_id': '',
                      'refseq_protein_id': '',
                      'openprot_id': '',
                      'uniprot_id': '',
                      'velia_id': velia_id_value,
                      'attrs': {}
                      }

    if include_table_name:
        orf_attributes = {'orf.' + k: v for k, v in orf_attributes.items()}

    logger.debug('Annotated ORF attributes: %s', orf_attributes)

    return orf_attributes, exon_indices, orf_cds_phases, orf_cds_frames


def validate_orf(orf_object: orf_classes.OrfBase, parent_transcript: dict, genome: Dict[str, str], accession_namespace: str) -> None:
    """Validate ORF.
    This function validates the ORF by comparing the internal ORF sequence with the extracted ORF sequence from the genome.
    Args:
        orf_object: The ORF object to validate.
        parent_transcript: The parent transcript of the ORF.
        genome: The genome from which the ORF is extracted.
        accession_namespace: The accession namespace for the genome sequences.

    Raises:
        AssertionError: If the internal ORF sequence does not match the extracted ORF sequence.
    """
    internal_orf_sequence = orf_object.nuc_sequence

    orf_attributes, exon_indices, orf_cds_phases, orf_cds_frames = populate_velia_db_fields_orf(
        orf_object, parent_transcript, include_table_name=True)
    assert len(exon_indices) == len(orf_cds_phases) == len(orf_cds_frames), f'ORF {orf_object}: length of exon indices {len(exon_indices)} does not equal length of CDS phases {len(orf_cds_phases)} or length of CDS frames {len(orf_cds_frames)}!'
    
    orf_attributes.update(parent_transcript)

    extracted_orf_sequence = sequence_tools.extract_orf_sequence_from_genome(
        orf_attributes, genome, source_coordinate_system='gff', accession_namespace=accession_namespace)

    assert len(
        extracted_orf_sequence) % 3 == 0, f'ORF {orf_object}: length {len(extracted_orf_sequence)} of extracted ORF sequence is not divisible by 3!'
    assert len(
        internal_orf_sequence) % 3 == 0, f'ORF {orf_object}: length {len(internal_orf_sequence)} of internal ORF sequence is not divisible by 3!'
    assert len(internal_orf_sequence) == len(
        extracted_orf_sequence), f'ORF {orf_object}: length {len(internal_orf_sequence)} of internal ORF sequence not equal to {len(extracted_orf_sequence)} of extracted ORF sequence!'
    assert internal_orf_sequence.upper() == extracted_orf_sequence.upper(
    ), f'ORF {orf_object}: internal ORF sequence {internal_orf_sequence} not equal to extracted ORF sequence {extracted_orf_sequence}!'


def generate_cds_orf_info(orf_attributes: Dict[str, Any],
                          parent_transcript: Dict[str, Any],
                          exon_indices: List[int],
                          cds_phases: List[int],
                          cds_frames: List[int],
                          exons_by_transcript: Dict[str, Dict[int, List[Dict[str, Any]]]],
                          cdss_by_exon: Dict[str, List[Dict[str, Any]]],
                          cds_orf_linkages: Set[Tuple[str, str]]) -> List[Dict[str, Any]]:
    """
    Generates CDS ORF information for a given ORF.

    This function iterates through exon indices, frames, and phases to generate rows of CDS ORF information. It logs warnings for exons with multiple definitions or no associated CDS entries.

    Args:
        orf_attributes: A dictionary containing attributes of the ORF.
        parent_transcript: A dictionary containing attributes of the parent transcript of the ORF.
        exon_indices: A list of exon indices relevant to the ORF.
        cds_phases: A list of phases for each CDS in the ORF.
        cds_frames: A list of reading frames for each CDS in the ORF.
        exons_by_transcript: A dictionary mapping transcript IDs to their exons.
        cdss_by_exon: A dictionary mapping exon IDs to their CDSs.

    Returns:
        A list of dictionaries, each representing a row of CDS ORF information to be inserted into a database or further processed.
    """
    logger.debug('Generating cds_orf entries for ORF %s in relation to transcript %s',
                 orf_attributes['orf.orf_idx_str'], parent_transcript['transcript.id'])
    transcript_exons = exons_by_transcript[parent_transcript['transcript.id']]

    cds_phases = misc.extract_numeric_vector(orf_attributes['orf.phases'])

    cds_orf_rows = []

    for exon_idx, cds_phase, cds_frame in zip(exon_indices, cds_phases, cds_frames):
        these_exons = transcript_exons[exon_idx + 1]
        if len(these_exons) > 1:
            logger.debug('Exon %s in transcript %s has multiple definitions!',
                         exon_idx, parent_transcript['transcript.id'])
        for this_exon in these_exons:
            if this_exon['exon.id'] in cdss_by_exon:
                for this_cds in cdss_by_exon[this_exon['exon.id']]:
                    if (this_cds['cds.id'], orf_attributes['orf.id']) in cds_orf_linkages: # This CDS-ORF linkage already exists, so we skip
                        continue

                    this_cds_orf_row = {'cds_orf.cds_id': this_cds['cds.id'], 'cds_orf.orf_id': orf_attributes['orf.id'],
                                        'cds_orf.orf_idx_str': orf_attributes['orf.orf_idx_str'],
                                        'cds_orf.cds_number': exon_idx, 
                                        'cds_orf.phase': cds_phase, 
                                        'cds_orf.reading_frame': cds_frame,
                                        'cds_orf.attrs': {}}

                    cds_orf_linkages.add((this_cds['cds.id'], orf_attributes['orf.id'])) # Add this CDS-ORF linkage to the set of linkages
                    cds_orf_rows.append(this_cds_orf_row)
            else:
                logger.debug(
                    'Exon %s has no associated CDS entries!', this_exon['exon.id'])

    return cds_orf_rows


def annotate_orfs(current_orf_id: int,
                  transcript_orfs: Dict[str, Set[orf_classes.OrfBase]],
                  transcripts_by_id: Dict[str, Dict[str, Union[str, int, float]]],
                  transcript_seq_strands_by_id: Dict[str, Tuple[str, str]],
                  orf_ids_by_idx_str: Dict[str, str],
                  genome: Dict[str, str],
                  exons_by_transcript: Dict[str, Dict[int, List[Dict[str, Union[str, int, float]]]]],
                  cdss_by_exon: Dict[str, List[Dict[str, Union[str, int, float]]]],
                  cds_orf_linkages: Set[Tuple[str, str]],
                  kozak_upstream_pssm: Bio.motifs.matrix.PositionSpecificScoringMatrix,
                  kozak_downstream_pssm: Bio.motifs.matrix.PositionSpecificScoringMatrix,
                  phase_style=constants.DEFAULT_PHASE_STYLE,
                  accession_namespace=constants.DEFAULT_ACCESSION_NAMESPACE,
                  skip_annotation: bool = False,
                  ) -> Tuple[
    int, List[Dict[str, Union[str, int, float]]], List[Dict[str, Union[str, int, float]]], List[
        Dict[str, Union[str, int, float]]], List[Dict[str, Union[str, int, float]]]]:
    """
    This function annotates Open Reading Frames (ORFs) in a single-threaded manner and returns comprehensive ORF-related data. It provides detailed ORF annotations, including basic ORF information, supplemental details, transcript-specific ORF data, and CDS-specific ORF data, alongside an incremented ORF ID.

    Args:
        current_orf_id (int): The identifier for the ORF currently being annotated.
        transcript_orfs (Dict[str, Set[orf_classes.OrfBase]]): Maps transcript identifiers to sets of ORF objects.
        transcripts_by_id (Dict[str, Dict[str, Union[str, int, float]]]): Maps transcript identifiers to their associated metadata, which may include a mix of strings, integers, and floats.
        transcript_seq_strands_by_id (Dict[str, Tuple[str, str]]): Maps transcript identifiers to tuples containing sequence data and strand orientation.
        genome (Dict[str, str]): A dictionary containing genome sequences indexed by their identifiers.
        exons_by_transcript (Dict[str, Dict[int, List[Dict[str, Union[str, int, float]]]]]): Maps transcript identifiers to their corresponding exons, which are represented as dictionaries of mixed data types.
        cdss_by_exon (Dict[str, List[Dict[str, Union[str, int, float]]]]): Maps exon identifiers to their corresponding CDS entries, represented as dictionaries of mixed data types.
        kozak_upstream_pssm (Bio.motifs.matrix.PositionSpecificScoringMatrix): A scoring matrix for analyzing the upstream region of the Kozak sequence for translation initiation.
        kozak_downstream_pssm (Bio.motifs.matrix.PositionSpecificScoringMatrix): A scoring matrix for analyzing the downstream region of the Kozak sequence for translation initiation.
        accession_namespace (str): The accession namespace for the genome sequences.

    Returns:
        Tuple[int, List[Dict[str, Union[str, int, float]]], List[Dict[str, Union[str, int, float]]], List[Dict[str, Union[str, int, float]]], List[Dict[str, Union[str, int, float]]]]: A tuple containing the updated ORF ID and four lists of dictionaries. Each list contains dictionaries with information about ORFs, supplemental ORFs, transcript ORFs, and CDS ORFs, respectively, featuring a mix of strings, integers, and floats.
    """
    logger.info('Annotating %s ORFs using a single process ...',
                sum((len(orfs) for orfs in transcript_orfs.values())))

    orf_table_chunk = []
    supplemental_orf_table_chunk = []
    transcript_orf_table_chunk = []
    cds_orf_table_chunk = []

    for transcript_id, orf_set in transcript_orfs.items():
        for orf_object in orf_set:
            parent_transcript = transcripts_by_id[transcript_id]
            parent_transcript_seq = transcript_seq_strands_by_id[transcript_id][0]
            logger.debug('Annotating %s ...', orf_object)

            validate_orf(orf_object=orf_object,
                         parent_transcript=parent_transcript,
                         genome=genome,
                         accession_namespace=accession_namespace)

            orf_row, exon_indices, cds_phases, cds_frames = populate_velia_db_fields_orf(orf_object, parent_transcript,
                                                                              str(current_orf_id),
                                                                              include_table_name=True,
                                                                              phase_style=phase_style)

            # Logic: We will often find the same ORF on multiple overlapping transcripts. These need to have a single
            # entry in the ORF table but multiple entries in the transcript_orf and cds_orf tables corresponding 
            # to the multiple transcripts / CDS sets on which the ORF is found.
            # So, the first ORF to have a given index string (which are identical for identical ORF coordinates)
            # will be assigned an ORF table entry and corresponding ORF ID number. This mapping between index strings
            # and ORF IDs is stored in orf_ids_by_idx_str. Subsequent ORFs with the same index
            # string will look up the ORF ID from orf_ids_by_idx_str and use that ORF ID for constructing the 
            # transcript_orf and cds_orf table entries.
            if orf_row['orf.orf_idx_str'] not in orf_ids_by_idx_str:
                orf_ids_by_idx_str[orf_row['orf.orf_idx_str']
                                   ] = orf_row['orf.id']
                logger.debug('ORF %s is the first at these coordinates. Assigned orf_id %s and adding to orf and supplemental_orf tables.', orf_row, orf_row['orf.id'])
                orf_table_chunk.append(orf_row)

                if not skip_annotation:
                    supplemental_orf_row = {'orf.id': current_orf_id,
                                            'orf.orf_idx_str': orf_row['orf.orf_idx_str'],
                                            'start_codon': orf_object.start_codon,
                                            'stop_codon': orf_object.stop_codon,
                                            'codon_length': orf_object.codon_length,
                                            'nuc_length': orf_object.nuc_length,
                                            'start_offset_from_transcript': orf_object.start_pos,
                                            }

                    supplemental_orf_row['kozak_score'] = motif_tools.compute_kozak_score(
                        orf_object=orf_object, parent_transcript_seq=parent_transcript_seq,
                        kozak_upstream_pssm=kozak_upstream_pssm,
                        kozak_downstream_pssm=kozak_downstream_pssm)

                current_orf_id += 1

            else: # This ORF is a duplicate, so we just need to add the transcript_orf and cds_orf entries
                logger.debug('ORF %s has already beenn called at these coordinates. Assigned existing orf_id %s. Will only add to transcript_orf and cds_orf tables.', orf_row, orf_ids_by_idx_str[orf_row['orf.orf_idx_str']])

                orf_row['orf.orf_id'] = orf_ids_by_idx_str[orf_row['orf.orf_idx_str']]

            if cdss_by_exon and exons_by_transcript:
                cds_orf_rows = generate_cds_orf_info(orf_attributes=orf_row,
                                                        parent_transcript=parent_transcript,
                                                        exon_indices=exon_indices,
                                                        cds_phases=cds_phases,
                                                        cds_frames=cds_frames,
                                                        exons_by_transcript=exons_by_transcript,
                                                        cdss_by_exon=cdss_by_exon,
                                                        cds_orf_linkages=cds_orf_linkages)
                cds_orf_table_chunk.extend(cds_orf_rows)               

            transcript_orf_row = {'transcript_orf.transcript_id': parent_transcript['transcript.id'],
                                  'transcript_orf.orf_id': orf_ids_by_idx_str[orf_row['orf.orf_idx_str']],
                                  'transcript_orf.orf_idx_str': orf_row['orf.orf_idx_str'],
                                  'transcript_orf.evidence_tag': 'big_prot',
                                  'transcript_orf.attrs': {}
                                  }
            transcript_orf_table_chunk.append(transcript_orf_row)

    logger.info('%s unique ORFs found in this chunk, %s total so far.',
                len(orf_table_chunk), len(orf_ids_by_idx_str))

    return current_orf_id, orf_table_chunk, supplemental_orf_table_chunk, transcript_orf_table_chunk, cds_orf_table_chunk


def convert_orf_row_to_object(orf_row_dict: Dict[str, Any], parent_transcript: Dict[str, Any], genome: Dict[str, str]) -> orf_classes.OrfBase:
    """
    Converts ORF row dictionary to ORF object.

    Args:
        orf_row_dict: The ORF row dictionary.
        parent_transcript: The parent transcript dictionary.
        genome: The genome dictionary.

    Returns:
        The ORF object.
    """
    logger.debug('Converting ORF %s to object ...',
                 orf_row_dict['orf.orf_idx_str'])

    orf_relative_start = coordinates.compute_orf_relative_start(
        orf_row_dict['orf.start'], parent_transcript)
    orf_sequence = sequence_tools.extract_orf_sequence_from_genome(
        orf_row_dict, genome).upper()
    assert len(orf_sequence) % 3 == 0, 'NT sequence for orf %s has length %s which is not evenly divisible by 3!' % (
        orf_row_dict['orf_idx_str'], len(orf_sequence))
    new_orf_object = orf_classes.OrfBase(parent_transcript['transcript.id'],
                                         start_pos=orf_relative_start,
                                         strand=orf_row_dict['orf.strand'],
                                         codons=sequence_tools.split_sequence_into_codons(
                                             orf_sequence)
                                         )
    return new_orf_object
