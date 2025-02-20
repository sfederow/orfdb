# OrfDB

OrfDB is a genomic annotation database that is designed to support open reading frames. 

## Overview

The OrfDB package provides a set of scripts and utilities to manage and load genomic data into a structured database (utilizing SQLAlchemy). This package contains the code to load various annotation sources into the ORF database, including GENCODE, RefSeq, CHESS, and custom transcriptomes. A utility to generate all possible ORFs from a set of transcripts (BigProt written by Dylan Skola) is also contained within this package. The database design is inspired by COBRAdb and the Ensembl genome database.

## Why OrfDB?

Standard genomic file formats (e.g. GFF) and databases (e.g. Ensembl) do not provide first-class support for ORFs. CDS and exons are well annotated and some support is provided for protein level entries but specific handling of ORFs is overlooked.  

## Installation

To install the OrfDB package, use the following command:

```
git clone https://github.com/veliatx/orfdb.git
pip install orfdb/
```

This will install the package and its dependencies, as well as set up the command line script `orfdb_load`.

## Data Download

Before running the database loading scripts, you'll need to download the required data files. Create a data directory and download the files as follows:

```bash
# Create data directory
mkdir -p orfdb/data
cd orfdb/data

# Download GENCODE annotation files
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.chr_patch_hapl_scaff.annotation.gff3.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.metadata.RefSeq.gz

# Download RefSeq files
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz

# Download CHESS files
wget https://github.com/chess-genome/chess/releases/download/v.3.1.3/chess3.1.3.GRCh38.gff.gz

# Download genome files
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_assembly_report.txt
```

Note: Replace [FTP_URL] with the appropriate FTP locations for each file. The specific versions of these files may be updated over time, so check the respective sources for the latest versions.

## Configuration

Before running the database loading scripts, ensure that the `settings.ini` file is properly configured. This file contains paths to the directories where the input files are located. Here is a brief overview of the configuration options:

- **DATABASE**: Configure the database connection settings, including host, port, user, password, and database name.
- **DATA**: Specify the filenames of the gff files downloaded above. These should all live in orfdb/data

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
