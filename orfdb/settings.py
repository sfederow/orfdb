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

# these are required
try:
    self.gencode_directory = Path(config['DATA']['gencode_directory'])
except NoOptionError:
    raise Exception('gencode_directory was not supplied in settings.ini')

try:
    self.gencode_version = config['DATA']['gencode_version']
except NoOptionError:
    raise Exception('gencode_version was not supplied in settings.ini')

try:
    self.genomes_directory = Path(config['DATA']['genomes_directory'])
except NoOptionError:
    raise Exception('genomes_directory was not supplied in settings.ini')

try:
    self.refseq_directory = Path(config['DATA']['refseq_directory'])
except NoOptionError:
    raise Exception('refseq_directory was not supplied in settings.ini')

try:
    self.uniprot_directory = Path(config['DATA']['uniprot_directory'])
except NoOptionError:
    raise Exception('uniprot_directory was not supplied in settings.ini')

try:
    self.chess_directory = Path(config['DATA']['chess_directory'])
except NoOptionError:
    raise Exception('chess_directory was not supplied in settings.ini')

try:
    self.openprot_directory = Path(config['DATA']['openprot_directory'])
except NoOptionError:
    raise Exception('openprot_directory was not supplied in settings.ini')

try:
    self.velia_directory = Path(config['DATA']['velia_directory'])
except NoOptionError:
    raise Exception('velia_directory was not supplied in settings.ini')

try:
    self.bigprot_directory = Path(config['DATA']['bigprot_directory'])
except NoOptionError:
    raise Exception('bigprot_directory was not supplied in settings.ini')

try:
    self.bigprot_version = get_highest_version_folder(self.bigprot_directory).name
except NoOptionError:
    raise Exception('bigprot_version was not supplied in settings.ini')
