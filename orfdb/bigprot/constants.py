"""
This module contains constants used in the ORF finding process. These include default values for minimum and maximum codon length, 
start and stop codons, and other parameters related to the ORF finding and annotation process.
"""
import os

DEFAULT_MIN_CODON_LENGTH = 15
DEFAULT_MAX_CODON_LENGTH = 1000
DEFAULT_START_CODONS = ('ATG', 'CTG', 'TTG', 'GTG')
DEFAULT_STOP_CODONS = ('TAA', 'TAG', 'TGA')
TRANSCRIPT_COORDINATE_SYSTEM = 'gff'
DEFAULT_PHASE_STYLE = 'gencode'
DEFAULT_ACCESSION_NAMESPACE = 'genbank'

DEFAULT_TRANSCRIPT_CHUNK_SIZE = 4096
DEFAULT_MP_CHUNK_SIZE = 128
DEFAULT_NUM_PROCESSES = os.cpu_count()

DEFAULT_PHYLO_CSF_GLOBS = [
    '/home/ubuntu/results/phylocsf_tracks/hg38_primates/final_tracks/PhyloCSFRaw*.bw']
DEFAULT_PHYLO_CSF_TRACKSET_NAMES = ['phylocsf_primates']

DEFAULT_GENBANK_FASTA_FPATH = '/home/ubuntu/data/hg38_genbank/GCA_000001405.28_GRCh38.p13_genomic.fna.gz'

DEFAULT_VELIADB_SETTINGS_FPATH = '/home/ubuntu/orfdb/settings.ini'

GFF_FIELDS = ['seqname', 'source', 'feature', 'start',
              'end', 'score', 'strand', 'frame', 'attribute']

DEFAULT_MOTIF_PSEUDOCOUNT = 0.5

DEFAULT_VALID_CHROMOSOMES = [str(i) for i in range(1, 24)] + ['X', 'Y']
DEFAULT_PROMOTER_UPSTREAM_WINDOW_SIZE = 30000
