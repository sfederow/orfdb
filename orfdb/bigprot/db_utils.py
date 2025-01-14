"""
This module provides utilities for interacting with the database.
"""

import configparser
import logging

import psycopg2

from . import constants

logger = logging.getLogger(__name__)
ORF_QUERY = '''SELECT
orf.aa_seq as "orf.aa_seq",
orf.assembly_id as "orf.assembly_id",
orf.attrs as "orf.attrs",
orf.benchling_id as "orf.benchling_id",
orf.block_sizes as "orf.block_sizes",
orf.chrom_starts as "orf.chrom_starts",
orf.end as "orf.end",
orf.ensembl_protein_id as "orf.ensembl_protein_id",
orf.id as "orf.id",
orf.nt_seq as "orf.nt_seq",
orf.openprot_id as "orf.openprot_id",
orf.orf_idx as "orf.orf_idx",
orf.orf_idx_str as "orf.orf_idx_str",
orf.phases as "orf.phases",
orf.refseq_protein_id as "orf.refseq_protein_id",
orf.secondary_orf_id as "orf.secondary_orf_id",
orf.start as "orf.start",
orf.strand as "orf.strand",
orf.uniprot_id as "orf.uniprot_id",
orf.velia_id as "orf.velia_id",
transcript_orf.transcript_id as "transcript_orf.transcript_id"
FROM orf
LEFT JOIN transcript_orf on orf.id = transcript_orf.orf_id
;'''

ORF_TRANSCRIPT_QUERY = '''SELECT DISTINCT ON (assembly.ucsc_style_name, orf.orf_idx_str)
assembly.assembly_unit as "assembly.assembly_unit",
assembly.assigned_molecule as "assembly.assigned_molecule",
assembly.assigned_molecule_location as "assembly.assigned_molecule_location",
assembly.attrs as "assembly.attrs",
assembly.genbank_accession as "assembly.genbank_accession",
assembly.genome_accession as "assembly.genome_accession",
assembly.id as "assembly.id",
assembly.refseq_accession as "assembly.refseq_accession",
assembly.sequence_length as "assembly.sequence_length",
assembly.sequence_role as "assembly.sequence_role",
assembly.ucsc_style_name as "assembly.ucsc_style_name",
orf.aa_seq as "orf.aa_seq",
orf.assembly_id as "orf.assembly_id",
orf.attrs as "orf.attrs",
orf.benchling_id as "orf.benchling_id",
orf.block_sizes as "orf.block_sizes",
orf.chrom_starts as "orf.chrom_starts",
orf.end as "orf.end",
orf.ensembl_protein_id as "orf.ensembl_protein_id",
orf.id as "orf.id",
orf.nt_seq as "orf.nt_seq",
orf.openprot_id as "orf.openprot_id",
orf.orf_idx as "orf.orf_idx",
orf.orf_idx_str as "orf.orf_idx_str",
orf.phases as "orf.phases",
orf.refseq_protein_id as "orf.refseq_protein_id",
orf.secondary_orf_id as "orf.secondary_orf_id",
orf.start as "orf.start",
orf.strand as "orf.strand",
orf.uniprot_id as "orf.uniprot_id",
orf.velia_id as "orf.velia_id",
transcript.assembly_id as "transcript.assembly_id",
transcript.attrs as "transcript.attrs",
transcript.block_sizes as "transcript.block_sizes",
transcript.chess_id as "transcript.chess_id",
transcript.chrom_starts as "transcript.chrom_starts",
transcript.end as "transcript.end",
transcript.ensembl_id as "transcript.ensembl_id",
transcript.gene_id as "transcript.gene_id",
transcript.id as "transcript.id",
transcript.refseq_id as "transcript.refseq_id",
transcript.start as "transcript.start",
transcript.strand as "transcript.strand",
transcript.support_level as "transcript.support_level",
transcript.transcript_idx as "transcript.transcript_idx",
transcript.transcript_idx_str as "transcript.transcript_idx_str",
transcript.velia_id as "transcript.velia_id",
transcript_orf.attrs as "transcript_orf.attrs",
transcript_orf.evidence_tag as "transcript_orf.evidence_tag",
transcript_orf.orf_id as "transcript_orf.orf_id",
transcript_orf.transcript_id as "transcript_orf.transcript_id"
FROM orf JOIN assembly on orf.assembly_id = assembly.id
JOIN transcript_orf on orf.id = transcript_orf.orf_id
JOIN transcript on transcript_orf.transcript_id = transcript.id
ORDER BY assembly.ucsc_style_name, orf_idx_str;'''

TRANSCRIPT_QUERY = '''SELECT DISTINCT ON (transcript.id)
transcript.id as "transcript.id", 
transcript.start as "transcript.start", 
transcript.end as "transcript.end", 
transcript.strand as "transcript.strand", 
transcript.chrom_starts as "transcript.chrom_starts", 
transcript.block_sizes as "transcript.block_sizes",
transcript.ensembl_id as "transcript.ensembl_id", 
transcript.refseq_id as "transcript.refseq_id",
transcript.transcript_type as "transcript.transcript_type",
assembly.id as "assembly.id",
assembly.genome_accession as "assembly.genome_accession",
assembly.ucsc_style_name as "assembly.ucsc_style_name",
assembly.genbank_accession as "assembly.genbank_accession",
assembly.sequence_length as "assembly.sequence_length"
FROM transcript
JOIN assembly ON assembly.id = transcript.assembly_id
;'''

ASSEMBLY_QUERY = '''SELECT
assembly.id as "assembly.id",
assembly.genbank_accession as "assembly.genbank_accession",
assembly.ucsc_style_name as "assembly.ucsc_style_name",
assembly.refseq_accession as "assembly.refseq_accession",
assembly.sequence_length as "assembly.sequence_length",
assembly.sequence_role as "assembly.sequence_role",
assembly.assigned_molecule as "assembly.assigned_molecule",
assembly.assigned_molecule_location as "assembly.assigned_molecule_location",
assembly.assembly_unit as "assembly.assembly_unit",
assembly.attrs as "assembly.attrs"
FROM assembly
ORDER BY assembly.id;'''

CDS_QUERY = '''SELECT
cds.id as "cds.id",
cds.ensembl_id as "cds.ensembl_id",
cds.refseq_id as "cds.refseq_id",
cds.chess_id as "cds.chess_id",
cds.openprot_id as "cds.openprot_id",
cds.velia_id as "cds.velia_id",
cds.ccds_id as "cds.ccds_id",
cds.ensembl_protein_id as "cds.ensembl_protein_id",
cds.refseq_protein_id as "cds.refseq_protein_id",
cds.attrs as "cds.attrs",
exon_cds.exon_id as "exon_cds.exon_id",
exon_cds.cds_id as "exon_cds.cds_id",
exon_cds.attrs as "exon_cds.attrs"
FROM cds JOIN exon_cds ON cds.id = exon_cds.cds_id;'''

EXON_QUERY = '''SELECT
exon.id as "exon.id",
exon.ensembl_id as "exon.ensembl_id",
exon.refseq_id as "exon.refseq_id",
exon.chess_id as "exon.chess_id",
exon.velia_id as "exon.velia_id",
exon.attrs as "exon.attrs",
transcript_exon.id as "transcript_exon.id",
transcript_exon.transcript_id as "transcript_exon.transcript_id",
transcript_exon.exon_id as "transcript_exon.exon_id",
transcript_exon.exon_number as "transcript_exon.exon_number",
transcript_exon.attrs as "transcript_exon.attrs"
FROM exon JOIN transcript_exon ON exon.id = transcript_exon.exon_id;'''


def filter_transcript_dict(transcripts_by_id: dict) -> dict:
    """
    Filters out transcripts with empty chrom_starts or block_sizes. So far this only seems to occur with 
    PAR-located genes on chrY, which are problematic for other reasons as well.

    Args:
        transcripts_by_id (dict): A dictionary of transcripts.

    Returns:
        dict: A filtered dictionary of transcripts.
    """
    filtered_transcripts = {tid: transcript for tid, transcript in transcripts_by_id.items() if
                            transcript['transcript.chrom_starts'] and transcript['transcript.block_sizes']}
    logger.info(
        'Removed %d transcripts, %d remain.', len(transcripts_by_id) - len(filtered_transcripts), len(filtered_transcripts))
    return filtered_transcripts


def describe_table(db_connection, table_name: str) -> list:
    """
    Describes a table in the database. Uses an empty query to get the column mames (since postgreSQL doesn't support "DESCRIBE").

    Args:
        db_connection: The database connection.
        table_name (str): The name of the table.

    Returns:
        list: A list of column names in the table.
    """
    cur = db_connection.cursor()
    cur.execute(f"""
    SELECT * FROM
    {table_name}
    WHERE false;
    """)

    return list([col.name for col in cur.description])


def fetch_one_as_dict(cursor, sort: bool = True) -> dict:
    """
    Fetches one row from the cursor as a dictionary.

    Args:
        cursor: The cursor to fetch from.
        sort (bool): Whether to sort the dictionary. Defaults to True.

    Returns:
        dict: The fetched row as a dictionary.
    """
    if sort:
        return {column.name: value for column, value in
                sorted(zip(cursor.description, cursor.fetchone()), key=lambda t: t[0])}
    return {column.name: value for column, value in zip(cursor.description, cursor.fetchone())}


def fetch_all_as_dict(cursor, sort: bool = True) -> list:
    """
    Fetches all rows from the cursor as a list of dictionaries.

    Args:
        cursor: The cursor to fetch from.
        sort (bool): Whether to sort the dictionaries. Defaults to True.

    Returns:
        list: The fetched rows as a list of dictionaries.
    """
    rows = []
    for row in cursor.fetchall():
        if sort:
            rows.append(
                {column.name: value for column, value in sorted(zip(cursor.description, row), key=lambda t: t[0])})
        else:
            rows.append({column.name: value for column,
                         value in zip(cursor.description, row)})
    return rows


def fetch_many_as_dict(cursor, size: int, sort: bool = True) -> list:
    """
    Fetches many rows from the cursor as a list of dictionaries.

    Args:
        cursor: The cursor to fetch from.
        size (int): The number of rows to fetch.
        sort (bool): Whether to sort the dictionaries. Defaults to True.

    Returns:
        list: The fetched rows as a list of dictionaries.
    """
    rows = []
    for row in cursor.fetchmany(size=size):
        if sort:
            rows.append(
                {column.name: value for column, value in sorted(zip(cursor.description, row), key=lambda t: t[0])})
        else:
            rows.append({column.name: value for column,
                         value in zip(cursor.description, row)})

    return rows


def get_config(settings_fpath: str = constants.DEFAULT_VELIADB_SETTINGS_FPATH) -> configparser.ConfigParser:
    """
    Gets the database configuration from a settings file.

    Args:
        settings_fpath (str): The path to the settings file. Defaults to constants.DEFAULT_VELIADB_SETTINGS_FPATH.

    Returns:
        ConfigParser: The database configuration.
    """
    with open(settings_fpath, 'rt', encoding='utf-8') as settings_file:
        db_config = configparser.ConfigParser()
        db_config.read_file(settings_file)
    return db_config


def connect(db_config: dict[str, dict[str, str]]) -> psycopg2.extensions.connection:
    return psycopg2.connect(database=db_config['DATABASE']['postgres_database'],
                            host=db_config['DATABASE']['postgres_host'],
                            user=db_config['DATABASE']['postgres_user'],
                            password=db_config['DATABASE']['postgres_password'],
                            port=db_config['DATABASE']['postgres_port'])


def query_engine_factory(database: str, host: str, user: str, password: str, port: str) -> callable:
    """
    Creates a query engine for executing queries on the database.

    Args:
        database (str): The name of the database.
        host (str): The host of the database.
        user (str): The user for the database.
        password (str): The password for the database.
        port (str): The port for the database.

    Returns:
        function: A function that executes queries on the database and returns the results.
    """

    def execute_query(query_string: str) -> list:
        conn = psycopg2.connect(database=database, host=host,
                                user=user, password=password,
                                port=port)
        cur = conn.cursor()
        cur.execute(query_string)
        results = fetch_all_as_dict(cur)
        conn.close()

        return results

    return execute_query


def make_veliadb_query_engine(veliadb_config_fpath: str = constants.DEFAULT_VELIADB_SETTINGS_FPATH) -> callable:
    """
    Creates a query engine for the Veliadb database.

    Args:
        veliadb_config_fpath (str): The path to the Veliadb configuration file. Defaults to constants.DEFAULT_VELIADB_SETTINGS_FPATH.

    Returns:
        function: A function that executes queries on the Veliadb database and returns the results.
    """
    veliadb_config = get_config(veliadb_config_fpath)
    return query_engine_factory(database=veliadb_config['DATABASE']['postgres_database'],
                                host=veliadb_config['DATABASE']['postgres_host'],
                                user=veliadb_config['DATABASE']['postgres_user'],
                                password=veliadb_config['DATABASE']['postgres_password'],
                                port=veliadb_config['DATABASE']['postgres_port'])


def get_assemblies(veliadb_config_fpath: str = constants.DEFAULT_VELIADB_SETTINGS_FPATH) -> list:
    """
    Gets all assemblies from the Veliadb database.

    Returns:
        list: A list of all assemblies in the Veliadb database.
    """
    logger.info("Retrieving assemblies from Velia DB ...")
    return make_veliadb_query_engine(veliadb_config_fpath)(ASSEMBLY_QUERY)


def get_orfs(veliadb_config_fpath: str = constants.DEFAULT_VELIADB_SETTINGS_FPATH) -> list:
    """
    Gets all ORFs from the Veliadb database.

    Returns:
        list: A list of all ORFs in the Veliadb database.
    """
    logger.info('Retrieving ORFs from Velia DB ...')
    return make_veliadb_query_engine(veliadb_config_fpath)(ORF_QUERY)


def get_transcripts(veliadb_config_fpath: str = constants.DEFAULT_VELIADB_SETTINGS_FPATH) -> list:
    """
    Gets all transcripts from the Veliadb database.

    Returns:
        list: A list of all transcripts in the Veliadb database.
    """
    logger.info('Retrieving transcripts from Veliadb ...')
    return make_veliadb_query_engine(veliadb_config_fpath)(TRANSCRIPT_QUERY)


def get_exons(veliadb_config_fpath: str = constants.DEFAULT_VELIADB_SETTINGS_FPATH) -> list:
    """
    Gets all exons from the Veliadb database.

    Returns:
        list: A list of all exons in the Veliadb database.
    """
    logger.info('Retrieving exons from Veliadb ...')

    return make_veliadb_query_engine(veliadb_config_fpath)(EXON_QUERY)



def get_cdss(veliadb_config_fpath: str = constants.DEFAULT_VELIADB_SETTINGS_FPATH) -> list:
    """
    Gets all CDSs from the Veliadb database.

    Returns:
        list: A list of all CDSs in the Veliadb database.
    """
    logger.info('Retrieving cdss from Veliadb ...')
    return make_veliadb_query_engine(veliadb_config_fpath)(CDS_QUERY)
