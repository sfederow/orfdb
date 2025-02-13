# -*- coding: utf-8 -*-

"""Retrive local user settings"""

from configparser import ConfigParser, NoOptionError
from pathlib import Path
from sys import modules
import os, re
from packaging.version import Version

self = modules[__name__]
settings_ini = (Path(*Path(os.path.abspath(__file__)).parts[:-2]) / 'settings.ini').resolve()

config_parser = ConfigParser()

def get_highest_version_folder(directory):
    """
    Finds the folder with the highest version number in the specified directory.
    
    Args:
        directory (str or Path): Path to the directory to scan.
    
    Returns:
        Path or None: The folder with the highest version, or None if no valid version folders are found.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError(f"{directory} is not a valid directory.")
    
    version_folders = {}
    version_pattern = re.compile(r"^v\d+(\.\d+)*$")  # Matches version-like patterns (e.g., 1.0.0)

    for folder in directory.iterdir():
        if folder.is_dir() and version_pattern.match(folder.name):
            try:
                version_folders[folder] = Version(folder.name)
            except ValueError:
                pass  # Ignore folders that don't follow proper versioning

    if version_folders:
        return max(version_folders, key=version_folders.get)
    return None


# overwrite defaults settings with settings from the file
if settings_ini.exists():
    config_parser.read(settings_ini)
    config = dict(config_parser)    
else:
    raise Exception('No settings files at path: %s' % settings_ini)

# set up the database connection string
self.db_connection_string = ('postgresql://%s:%s@%s/%s' %
                             (config['DATABASE']['postgres_user'],
                              config['DATABASE']['postgres_password'],
                              config['DATABASE']['postgres_host'],
                              config['DATABASE']['postgres_database']
                            ))

self.data_dir = Path(settings_ini).parent.joinpath('data')

try:
    self.gencode_gff = Path(config['DATA']['gencode_gff'])
except NoOptionError:
    raise Exception('gencode_gff was not supplied in settings.ini')

try:
    self.gencode_refseq = Path(config['DATA']['gencode_refseq'])
except NoOptionError:
    raise Exception('gencode_refseq was not supplied in settings.ini')

try:
    self.genome_assembly = Path(config['DATA']['genome_assembly'])
except NoOptionError:
    raise Exception('genomes_assembly was not supplied in settings.ini')

try:
    self.genome = Path(config['DATA']['genome'])
except NoOptionError:
    raise Exception('genome was not supplied in settings.ini')

try:
    self.refseq_gff = Path(config['DATA']['refseq_gff'])
except NoOptionError:
    raise Exception('refseq_gff was not supplied in settings.ini')

try:
    self.chess_gff = Path(config['DATA']['chess_gff'])
except NoOptionError:
    raise Exception('chess_gff was not supplied in settings.ini')

try:
    self.bigprot_directory = Path(settings_ini).parent.joinpath('data', 'bigprot')
except NoOptionError:
    raise Exception('bigprot_directory was not supplied in settings.ini')

try:
    self.bigprot_version = get_highest_version_folder(self.bigprot_directory)
except NoOptionError:
    raise Exception('bigprot_version was not supplied in settings.ini')
