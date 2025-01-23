# OrfDB

OrfDB is a genomic annotation database that is designed to support open reading frames. 

## Overview

The OrfDB package provides a set of scripts and utilities to manage and load genomic data into a structured database. This package contains the code to load various annotation sources into the ORF database, including GENCODE, RefSeq, CHESS, and custom transcriptomes. A utility to generate all possible ORFs from a set of transcripts (BigProt written by Dylan Skola) is also contained within this package. The database design is inspired by COBRAdb and the Ensembl genome database.

## Why OrfDB?

Standard genomic file formats (e.g. GFF) and databases (e.g. Ensembl) do not provide first-class support for ORFs. CDS and exons are well annotated and some support is provided for protein level entries but specific handling of ORFs is overlooked.  

## Installation

To install the OrfDB package, use the following command:

```
git clone https://github.com/veliatx/orfdb.git
pip install orfdb/
```

This will install the package and its dependencies, as well as set up the command line script `orfdb_load`.

## Configuration

Before running the database loading scripts, ensure that the `settings.ini` file is properly configured. This file contains paths to the directories where the input files are located. Here is a brief overview of the configuration options:

- **DATABASE**: Configure the database connection settings, including host, port, user, password, and database name.
- **DATA**: Specify the directories containing the various annotation files, such as GENCODE, RefSeq, CHESS, and others.

## Running the Database Loading Script

The main entry point for loading data into the ORF database is the `load_db` function, which can be executed via the command line using the installed script:

```
orfdb_load
```

### Command Line Options

- `--drop-all`: Use this flag to empty the database and reload all data. This is useful for resetting the database to a clean state before loading new data.

Example usage:

```
orfdb_load --drop-all
```

This command will drop all existing tables in the database and recreate them before loading the data.

## Logging

The script will generate a log file with the name format `YYYYMMDD_HHMMSS_orfdb_load_db.log`, which contains detailed information about the loading process, including any errors encountered.

## Contributing

Contributions to OrfDB are welcome. Please ensure that any changes are well-documented and tested.

## License

OrfDB is licensed under the MIT License. See the LICENSE file for more details.

## Contact

For any questions or issues, please contact the author:

- **Author**: Stephen Federowicz, Dylan Skola
- **Email**: sfederow@gmail.com, dylan@phageghost.net
