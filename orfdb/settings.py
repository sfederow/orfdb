# -*- coding: utf-8 -*-

"""Retrive local user settings"""

from configparser import ConfigParser, NoOptionError
from pathlib import Path
from sys import modules
import os

self = modules[__name__]
settings_ini = (Path(*Path(os.path.abspath(__file__)).parts[:-2]) / 'settings.ini').resolve()

config_parser = ConfigParser()

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