"""
This script is used to find Open Reading Frames (ORFs) in a set of sequences. It takes as input a FASTA file containing genomic sequences, 
and a settings.ini file containing the credentials to access the veliaDB. The script also allows for customization of codon settings, 
including start and stop codons, and minimum and maximum codon length. Advanced settings such as the glob pattern for files containing 
phylocsf bigwig tracks for 3 frames on both strands, and the number of simultaneous processes to run can also be specified. 
The output is logged in a specified output path.
"""

import argparse
import gzip
import logging
import os, re
import sys
from typing import List

import Bio.SeqIO
from tqdm import tqdm

from orfdb.bigprot import bigwig_tracks
from orfdb.bigprot import constants
from orfdb.bigprot import db_utils
from orfdb.bigprot import gwas_tools
from orfdb.bigprot import gff_tools
from orfdb.bigprot import orf_calling
from orfdb import settings

logger = logging.Logger('placeholder')


def setup_logger(log_fpath: str, verbosity: int = logging.INFO) -> None:
    """
    Set up the logger.

    Args:
        log_fpath (str): The path to the log file.
        verbosity (int, optional): The level of logging. Defaults to logging.INFO.
    """
    global logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(verbosity)

    # create file handler
    fh = logging.FileHandler(log_fpath, mode='w')
    fh.setLevel(verbosity)

    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    # add the handlers to logger
    logger.addHandler(ch)
    logger.addHandler(fh)

    logger.info('Logger set up at level %s and writing to %s.', 
        verbosity, log_fpath)


def setup_argparser() -> argparse.Namespace:
    """
    Set up and execute the argument parser.

    Returns:
        argparse.Namespace: The parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description='Find ORFs in a set of sequences.')

    # Group: Basic settings
    basic_group = parser.add_argument_group('Basic settings')

    basic_group_exclusive = basic_group.add_mutually_exclusive_group(
        required=True)
    basic_group_exclusive.add_argument('-s', '--db-settings', type=str,
                                       help='Path to the settings.ini file that contains the credentials to access veliaDB for obtaining transcript information. Mutally exlusive with --gtf-file.',)
    basic_group_exclusive.add_argument('-gf', '--gtf-file', type=str,
                                       help='Path to a GTF file from which to obtain transcript information. Mutally exclusive with --db-settings.')
    basic_group.add_argument('-o', '--output', type=str, required=True,
                             help='Output path.')
    basic_group.add_argument('-v', '--verbose', action='store_true',
                             help='Set verbose logging on.')
    basic_group.add_argument('-d', '--dataset-name', type=str, required=True,
                             help='Name of the dataset.')
    basic_group.add_argument('-f', '--fasta-file', type=str, required=True,
                             help='Path to a FASTA file containing genomic sequences.')

    # Group: Codon settings
    codon_group = parser.add_argument_group('Codon settings')
    codon_group.add_argument('-sc', '--start-codons', type=str, default=constants.DEFAULT_START_CODONS,
                             help=f'Start codons (optional, default value {",".join(constants.DEFAULT_START_CODONS)}).')
    codon_group.add_argument('-stc', '--stop-codons', type=str, default=constants.DEFAULT_STOP_CODONS,
                             help=f'Stop codons (optional, default value {",".join(constants.DEFAULT_STOP_CODONS)}).')
    codon_group.add_argument('-min', '--min-codon-length', type=int, default=constants.DEFAULT_MIN_CODON_LENGTH,
                             help=f'Minimum codon length for called ORFs (optional, default value {constants.DEFAULT_MIN_CODON_LENGTH}).')
    codon_group.add_argument('-max', '--max-codon-length', type=int, default=constants.DEFAULT_MAX_CODON_LENGTH,
                             help=f'Maximum codon length for called ORFs (optional, default value {constants.DEFAULT_MAX_CODON_LENGTH}).')

    annotation_group = parser.add_argument_group('Annotation settings')
    annotation_group.add_argument('-sa', '--skip-annotation', action='store_true',
                                  help='Skip annotation of ORFs, just call them.')
    annotation_group.add_argument('-g', '--gwas-associations-fpath', type=str, required=False,
                                  help='Path to a TSV file containing SNP associations (in the format of NHGRI-EBI)')
    annotation_group.add_argument('-p', '--phylocsf-globs', type=str, nargs='+', required=False,
                                  default=[],
                                  help='Glob patterns for files containing phylocsf bigwig tracks for 3 frames on both strands.')
    annotation_group.add_argument('-pn', '--phylocsf-trackset-names', type=str, nargs='+', required=False,
                                  default=[],
                                  help='Names for the phylocsf tracksets. Will be prepended to the names of the summary stats.')
    # Group: Advanced settings
    advanced_group = parser.add_argument_group('Advanced settings')

    advanced_group.add_argument('-np', '--num-processes', type=int, default=constants.DEFAULT_NUM_PROCESSES,
                                help='Number of simultaneous processes to run (optional, default to the max available).')
    advanced_group.add_argument('-tc', '--transcript-chunk-size', type=int,
                                default=constants.DEFAULT_TRANSCRIPT_CHUNK_SIZE,
                                help=f'The size of transcript chunks for processing (optional, default value {constants.DEFAULT_TRANSCRIPT_CHUNK_SIZE}).')
    advanced_group.add_argument('-ps', '--phase-style', type=str,
                                default=constants.DEFAULT_PHASE_STYLE,
                                help='The phase style to use for CDSs. Can be "gencode" or "ensembl".')
    advanced_group.add_argument('-an', '--accession-namespace', type=str, default=constants.DEFAULT_ACCESSION_NAMESPACE,
                                help=f'The accession namespace for the genome sequences. Typical values are "ensembl", "genbank" (optional, default value "{constants.DEFAULT_ACCESSION_NAMESPACE}").')
    advanced_group.add_argument('-mc', '--max-chunks', type=int, default=-1,
                                help='Maximum number of transcript chunks to process. Useful to launch small runs for debugging purposes. If not specified, all transcript chunks will be processed.')
    args = parser.parse_args()

    return args


def perform_analysis(output: str,
                     dataset_name: str,
                     genome_fasta_fpath: str,
                     db_settings_fpath: str,
                     skip_annotation: bool,
                     min_codon_length: int,
                     max_codon_length: int,
                     transcript_chunk_size: int,
                     max_chunks: int,
                     num_processes: int,
                     verbose: bool,
                     accession_namespace: str,
                     phase_style=constants.DEFAULT_PHASE_STYLE,
                     version: str = None) -> None:
    """
    Conducts the analysis by loading and processing the specified datasets. This includes identifying and optionally annotating Open Reading Frames (ORFs), handling genomic sequences, and analyzing SNP associations along with PhyloCSF bigwig tracks for evolutionary conservation insights.

    Args:
        output (str): The directory path where output files will be saved.
        dataset_name (str): The name identifying the dataset being analyzed.
        genome_fasta_fpath (str): The file path to the FASTA file containing genomic sequences.
        db_settings_fpath (str): The file path to the database settings file (settings.ini) with configuration details.
        min_codon_length (int): The minimum codon length for identifying ORFs.
        max_codon_length (int): The maximum codon length for identifying ORFs.
        transcript_chunk_size (int): The size of chunks into which transcripts are divided for processing.
        max_chunks (int): The maximum number of transcript chunks to process, useful for debugging or limited runs.
        num_processes (int): The number of processes to run concurrently during the analysis.
        verbose (bool): If true, enables verbose logging for more detailed execution logs.
        skip_annotation (bool): If true, skips the annotation of ORFs, focusing only on their identification.
        phase_style (str): The phase style to be used for CDSs, indicating the coding sequence phase convention ('gencode' or 'ensembl').
        accession_namespace (str): The namespace for the accession numbers of the genome sequences.
    """
    os.makedirs(output, exist_ok=True)
    output_prefix = os.path.join(output,
                                 f'orfset_{dataset_name}_minlen_{min_codon_length}_maxlen_{max_codon_length}')
    log_file_fpath = f'{output_prefix}.find_orfs.log'

    if verbose:
        setup_logger(log_fpath=log_file_fpath,
                     verbosity=logging.DEBUG)
    else:
        setup_logger(log_fpath=log_file_fpath,
                     verbosity=logging.INFO)


    #version = settings.bigprot_version.name
    #version = bump_version(version)

    # Log the function arguments
    logger.info(f"Calling and annotating ORFs using big_prot version {version} with the following arguments: "
                f"output={output}, dataset_name={dataset_name}, genome_fasta_fpath={genome_fasta_fpath}, "
                f"settings={db_settings_fpath}, "
                f"min_codon_length={min_codon_length}, max_codon_length={max_codon_length}, "
                f"transcript_chunk_size={transcript_chunk_size}, max_chunks={max_chunks}, num_processes={num_processes}, "
                f"verbose={verbose}, skip_annotation={skip_annotation}, phase_style={phase_style}, "
                f"accession_namespace={accession_namespace}.")


    # Load the genome as a dictionary of Bio.SeqRecord objects keyed by accession. ToDo: Standardize on either SeqRecords or strings for sequences!
    logger.info(f'Loading genome sequences from {genome_fasta_fpath} ...')
    with gzip.open(genome_fasta_fpath, 'rt') as fasta_file:
        genome = {seq_record.id: str(seq_record.seq).upper() for seq_record in tqdm(
            Bio.SeqIO.parse(fasta_file, format='fasta'))}
    logger.info(f'Loaded {len(genome)} sequences.')

    exons_by_transcript = {}
    cdss_by_exon = {}

    if db_settings_fpath is None:
        logger.info('No DB settings file provided, skipping DB connection.')
        
    else:
        # Load database tables
        logger.info('Connecting to DB defined in settings file %s',
                    db_settings_fpath)
        execute_veliadb_query = db_utils.make_veliadb_query_engine(
            db_settings_fpath)
        logger.info('Querying transcript table ...')
        transcripts = execute_veliadb_query(db_utils.TRANSCRIPT_QUERY)
        logger.info('Received %d transcripts.', len(transcripts))
        # Turn into dictionary and remove PAR genes
        transcripts_by_id = db_utils.filter_transcript_dict(
            {str(transcript['transcript.id']): transcript for transcript in tqdm(transcripts)})
        del transcripts

        logger.info('Querying exons ...')

        exon_count = 0
        for exon in tqdm(execute_veliadb_query(db_utils.EXON_QUERY)):
            exon_num = int(exon['transcript_exon.exon_number'])
            if exon['transcript_exon.transcript_id'] not in exons_by_transcript:
                exons_by_transcript[exon['transcript_exon.transcript_id']] = {}
            if exon_num not in exons_by_transcript[exon['transcript_exon.transcript_id']]:
                exons_by_transcript[exon['transcript_exon.transcript_id']][exon_num] = [
                ]
            exons_by_transcript[exon['transcript_exon.transcript_id']][exon_num].append(
                exon)
            exon_count += 1
        logger.info('Received %d exons.', exon_count)

        logger.info('Querying CDSs ...')
        cds_count = 0
        for cds in tqdm(execute_veliadb_query(db_utils.CDS_QUERY)):
            if cds['exon_cds.exon_id'] not in cdss_by_exon:
                cdss_by_exon[cds['exon_cds.exon_id']] = []
            cdss_by_exon[cds['exon_cds.exon_id']].append(cds)
            cds_count += 1
        logger.info('Received %d CDSs.', cds_count)

    orf_calling.call_and_annotate_orfs(transcripts_by_id=transcripts_by_id,
                                       genome=genome,
                                       output_prefix=output_prefix,
                                       exons_by_transcript=exons_by_transcript,
                                       cdss_by_exon=cdss_by_exon,
                                       skip_annotation=skip_annotation,
                                       min_codon_length=min_codon_length,
                                       max_codon_length=max_codon_length,
                                       phase_style=phase_style,
                                       accession_namespace=accession_namespace,
                                       transcript_chunk_size=transcript_chunk_size,
                                       max_chunks=max_chunks,
                                       num_processes=num_processes)


def main() -> None:
    """
    Main function to parse command line arguments, and initiate the analysis.
    """
    args = setup_argparser()

    print(args)

    perform_analysis(output=args.output, verbose=args.verbose, dataset_name=args.dataset_name,
                     genome_fasta_fpath=args.fasta_file,
                     db_settings_fpath=args.db_settings,
                     gtf_file_fpath=args.gtf_file,
                     skip_annotation=args.skip_annotation,
                     min_codon_length=args.min_codon_length,
                     max_codon_length=args.max_codon_length,
                     phase_style=args.phase_style,
                     accession_namespace=args.accession_namespace,
                     transcript_chunk_size=args.transcript_chunk_size,
                     max_chunks=args.max_chunks,
                     num_processes=args.num_processes)


if __name__ == '__main__':
    sys.exit(main())
