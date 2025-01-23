# OrfDB

OrfDB is a comprehensive resource for internal small open reading frame (sORF) annotation and experimental data. This package is designed to load various annotation sources into the ORF database, including GENCODE, RefSeq, CHESS, OpenProt, and custom ORF annotations.

## Overview

The OrfDB package provides a set of scripts and utilities to manage and load genomic data into a structured database. It supports multiple data sources and formats, ensuring that the database is populated with the most relevant and up-to-date information.

## Installation

To install the OrfDB package, use the following command:

bash
pip install .

This will install the package and its dependencies, as well as set up the command line script `orfdb_load`.

## Configuration

Before running the database loading scripts, ensure that the `settings.ini` file is properly configured. This file contains paths to the directories where the input files are located. Here is a brief overview of the configuration options:

- **DATABASE**: Configure the database connection settings, including host, port, user, password, and database name.
- **DATA**: Specify the directories containing the various annotation files, such as GENCODE, RefSeq, CHESS, and others.

## Running the Database Loading Script

The main entry point for loading data into the ORF database is the `load_db` function, which can be executed via the command line using the installed script:

bash
orfdb_load


### Command Line Options

- `--drop-all`: Use this flag to empty the database and reload all data. This is useful for resetting the database to a clean state before loading new data.

Example usage:

bash
orfdb_load --drop-all


This command will drop all existing tables in the database and recreate them before loading the data.

## Logging

The script will generate a log file with the name format `YYYYMMDD_HHMMSS_orfdb_load_db.log`, which contains detailed information about the loading process, including any errors encountered.

## Contributing

Contributions to OrfDB are welcome. Please ensure that any changes are well-documented and tested.

## License

OrfDB is licensed under the MIT License. See the LICENSE file for more details.

## Contact

For any questions or issues, please contact the author:

- **Author**: Stephen Federowicz
- **Email**: steve@veliatx.com
