"""Module to implement ORM for the Velia Database"""

from veliadb.settings import db_connection_string
from os import system

from sqlalchemy.orm import sessionmaker, relationship, aliased
from sqlalchemy.orm.session import Session as _SA_Session
from sqlalchemy import Table, MetaData, create_engine, Column, Integer, \
    String, Boolean, ForeignKey, and_, or_, not_, distinct, select, Identity, JSON
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import insert
from types import MethodType


engine = create_engine(db_connection_string)
Base = declarative_base(bind=engine)
metadata = MetaData(bind=engine)


class Assembly(Base):
    """Represents a genomic assembly in the database.

    Stores information about genome assemblies including accession numbers,
    sequence details, and assembly metadata.

    Attributes:
        id (int): Primary key
        genbank_accession (str): GenBank accession number
        refseq_accession (str): RefSeq accession number
        ucsc_style_name (str): UCSC genome browser style name
        sequence_role (str): Role of the sequence
        assembly_unit (str): Assembly unit identifier
        assigned_molecule (str): Assigned molecule name
        assigned_molecule_location (str): Location on assigned molecule
        sequence_length (int): Length of the sequence
        genome_accession (str): Genome accession number
        attrs (JSON): Additional attributes stored as JSON
    """
    __tablename__ = 'assembly'

    id = Column(Integer, primary_key=True, autoincrement=True)
    genbank_accession = Column(String(200))
    refseq_accession = Column(String(100))
    ucsc_style_name = Column(String(50))
    sequence_role = Column(String(200))
    assembly_unit = Column(String(200))
    assigned_molecule = Column(String(200))
    assigned_molecule_location = Column(String(200))
    sequence_length = Column(Integer, nullable=False)
    genome_accession = Column(String(200))
    attrs = Column(JSON)

    __table_args__ = (UniqueConstraint('genbank_accession', 'refseq_accession'),{})

    def __repr__(self):
        return f"Assembly (#{self.id}): {self.genome_accession}, {self.ucsc_style_name},"\
               f"{self.genbank_accession}, {self.sequence_length}"

    def __init__(self, genbank_accession, refseq_accession, ucsc_style_name,
                 sequence_role, assembly_unit, assigned_molecule, 
                 assigned_molecule_location, sequence_length, genome_accession):
        self.genbank_accession = genbank_accession
        self.refseq_accession = refseq_accession
        self.ucsc_style_name = ucsc_style_name
        self.sequence_role = sequence_role
        self.assembly_unit = assembly_unit
        self.assigned_molecule = assigned_molecule
        self.assigned_molecule_location = assigned_molecule_location
        self.sequence_length = sequence_length
        self.genome_accession = genome_accession


class SequenceRegion(Base):
    """Base class for genomic sequence regions.

    Represents a region of sequence with coordinates and strand information.
    Parent class for Gene, Exon, UtrThree, UtrFive, and Cds.

    Attributes:
        id (int): Primary key
        assembly_id (int): Foreign key to Assembly
        start (int): Start position
        end (int): End position
        strand (str): Strand ('+' or '-')
        type (str): Type of sequence region
    """
    __tablename__ = 'sequence_region'

    id = Column(Integer, primary_key=True, autoincrement=True)
    assembly_id = Column(Integer, ForeignKey('assembly.id'))
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    strand = Column(String(1), nullable=False)
    type = Column(String(30))

    __table_args__ = (UniqueConstraint('start', 'end', 'strand', 'assembly_id', 'type'),{})

    __mapper_args__ = {'polymorphic_identity': 'sequence_region',
                       'polymorphic_on': type
                      }

    def __repr__(self):
        return f"SequenceRegion: {self.start}-{self.end} ({self.strand})"

    def __repr__dict__(self):
        return {"id":self.id,"start":self.start,"end":self.end,"strand":self.strand}

    def __init__(self, start, end, strand, assembly_id, type):
        self.start = start
        self.end = end
        self.strand = strand
        self.assembly_id = assembly_id
        self.type = type


class Gene(SequenceRegion):
    """Represents a gene annotation.

    Stores gene information including various database identifiers and
    relationships to transcripts.

    Attributes:
        id (int): Primary key, inherited from SequenceRegion
        hgnc_id (str): HGNC identifier
        hgnc_name (str): HGNC gene name
        ensembl_id (str): Ensembl gene ID
        refseq_id (str): RefSeq gene ID
        chess_id (str): CHESS database ID
        velia_id (str): Velia database ID
        gene_type (str): Type of gene
        long_name (str): Full gene name
        attrs (JSON): Additional attributes
        transcripts (relationship): Related Transcript objects
    """
    __tablename__ = 'gene'

    id = Column(Integer, ForeignKey('sequence_region.id'), primary_key=True)
    hgnc_id = Column(String(30))
    hgnc_name = Column(String(30))
    ensembl_id = Column(String(30))
    refseq_id = Column(String(30))
    chess_id = Column(String(30))
    velia_id = Column(String(30))
    gene_type = Column(String(200))
    long_name = Column(String(100))
    attrs = Column(JSON)

    transcripts = relationship('Transcript', primaryjoin='Gene.id==Transcript.gene_id')

    __mapper_args__ = { 'polymorphic_identity': 'gene' }

    def __repr__(self):
        return f"Gene: ({self.hgnc_id}, {self.ensembl_id}, {self.hgnc_name})"\
               f"{self.start}-{self.end} ({self.strand})"

    def __init__(self, start, end, strand, assembly_id, hgnc_id, hgnc_name, 
                 ensembl_id, refseq_id, chess_id, velia_id, gene_type, long_name=None, attrs={}):
        
        super(Gene, self).__init__(start, end, strand, assembly_id, type='gene')
        self.hgnc_id = hgnc_id
        self.hgnc_name = hgnc_name
        self.ensembl_id = ensembl_id
        self.refseq_id = refseq_id
        self.chess_id = chess_id
        self.velia_id = velia_id
        self.gene_type = gene_type
        self.long_name = long_name
        self.attrs = attrs


class Transcript(Base):
    """Represents a transcript annotation.

    Stores transcript information including coordinates, identifiers,
    and relationships to genes and exons.

    Attributes:
        id (int): Primary key
        start (int): Start position
        end (int): End position
        strand (str): Strand ('+' or '-')
        assembly_id (int): Foreign key to Assembly
        block_sizes (str): Sizes of transcript blocks
        chrom_starts (str): Chromosome start positions
        transcript_idx (str): Unique transcript index
        transcript_idx_str (str): Human-readable transcript index
        gene_id (int): Foreign key to Gene
        ensembl_id (str): Ensembl transcript ID
        refseq_id (str): RefSeq transcript ID
        chess_id (str): CHESS database ID
        velia_id (str): Velia database ID
        support_level (str): Transcript support level
        transcript_type (str): Type of transcript
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'transcript'

    id = Column(Integer, primary_key=True, autoincrement=True)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    strand = Column(String(1), nullable=False)
    assembly_id = Column(Integer, ForeignKey('assembly.id'))
    block_sizes = Column(String)
    chrom_starts = Column(String)
    transcript_idx = Column(String(130))
    transcript_idx_str = Column(String)

    #utr_five_id = Column(Integer, ForeignKey('utr_five.id'))
    #utr_three_id = Column(Integer, ForeignKey('utr_three.id'))

    gene_id = Column(Integer, ForeignKey('gene.id'))
    ensembl_id = Column(String(30))
    refseq_id = Column(String(30))
    chess_id = Column(String(30))
    velia_id = Column(String(30))
    support_level = Column(String)
    transcript_type = Column(String)
    attrs = Column(JSON)

    gene = relationship('Gene', primaryjoin='Transcript.gene_id==Gene.id', viewonly=True)
    exons = relationship('Exon', secondary='transcript_exon', viewonly=True)

    __table_args__ = (UniqueConstraint('transcript_idx'), {})

    def __repr__(self):
        return f"Transcript: ({self.id}, {self.ensembl_id}, {self.refseq_id})"\
               f"{self.strand}:{self.start}-{self.end}"

    def __init__(self, start, end, strand, assembly_id, block_sizes, chrom_starts, 
                 transcript_idx, transcript_idx_str, gene_id, ensembl_id, 
                 refseq_id, chess_id, velia_id, support_level, transcript_type, attrs={}):
        
        self.start = start
        self.end = end
        self.strand = strand
        self.assembly_id = assembly_id
        self.block_sizes = block_sizes
        self.chrom_starts = chrom_starts
        self.transcript_idx = transcript_idx
        self.transcript_idx_str = transcript_idx_str

        self.gene_id = gene_id
        self.ensembl_id = ensembl_id
        self.refseq_id = refseq_id
        self.chess_id = chess_id
        self.velia_id = velia_id
        self.support_level = support_level
        self.transcript_type = transcript_type
        self.attrs = attrs


class Exon(SequenceRegion):
    """Represents an exon annotation.

    Stores exon information and relationships to transcripts.

    Attributes:
        id (int): Primary key, inherited from SequenceRegion
        ensembl_id (str): Ensembl exon ID
        refseq_id (str): RefSeq exon ID
        chess_id (str): CHESS database ID
        velia_id (str): Velia database ID
        attrs (JSON): Additional attributes
        transcripts (relationship): Related Transcript objects
    """
    __tablename__ = 'exon'

    id = Column(Integer, ForeignKey('sequence_region.id'), primary_key=True)
    ensembl_id  = Column(String(30))
    refseq_id = Column(String(30))
    chess_id = Column(String(30))
    velia_id = Column(String(30))
    attrs = Column(JSON)

    transcripts = relationship('Transcript', secondary='transcript_exon', viewonly=True)

    __mapper_args__ = { 'polymorphic_identity': 'exon' }

    def __repr__(self):
        return f"Exon: ({self.id}, {self.ensembl_id}) {self.start}-{self.end}"

    def __init__(self, start, end, strand, assembly_id, ensembl_id, refseq_id, 
                chess_id, velia_id, attrs={}):
       
        super(Exon, self).__init__(start, end, strand, assembly_id, 'exon')
        self.ensembl_id = ensembl_id
        self.refseq_id = refseq_id
        self.chess_id = chess_id
        self.velia_id = velia_id
        self.velia_id = velia_id
        self.attrs = attrs


class UtrThree(SequenceRegion):
    """Represents a 3' UTR annotation.

    Stores information about 3' untranslated regions.

    Attributes:
        id (int): Primary key, inherited from SequenceRegion
        ensembl_id (str): Ensembl UTR ID
        refseq_id (str): RefSeq UTR ID
        chess_id (str): CHESS database ID
        velia_id (str): Velia database ID
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'utr_three'

    id = Column(Integer, ForeignKey('sequence_region.id'), primary_key=True)

    ensembl_id  = Column(String(30))
    refseq_id = Column(String(30))
    chess_id = Column(String(30))
    velia_id = Column(String(30))
    attrs = Column(JSON)

    __mapper_args__ = { 'polymorphic_identity': 'utr_three' }

    def __repr__(self):
        return f"3' UTR: ({self.id}, {self.ensembl_id}) {self.start}-{self.end}"

    def __init__(self, start, end, strand, assembly_id, ensembl_id, refseq_id, 
                chess_id, velia_id, attrs={}):
       
        super(UtrThree, self).__init__(start, end, strand, assembly_id, 'utr_three')
        self.ensembl_id = ensembl_id
        self.refseq_id = refseq_id
        self.chess_id = chess_id
        self.velia_id = velia_id
        self.attrs = attrs


class UtrFive(SequenceRegion):
    """Represents a 5' UTR annotation.

    Stores information about 5' untranslated regions.

    Attributes:
        id (int): Primary key, inherited from SequenceRegion
        ensembl_id (str): Ensembl UTR ID
        refseq_id (str): RefSeq UTR ID
        chess_id (str): CHESS database ID
        velia_id (str): Velia database ID
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'utr_five'

    id = Column(Integer, ForeignKey('sequence_region.id'), primary_key=True)

    ensembl_id  = Column(String(30))
    refseq_id = Column(String(30))
    chess_id = Column(String(30))
    velia_id = Column(String(30))
    attrs = Column(JSON)

    __mapper_args__ = { 'polymorphic_identity': 'utr_five' }

    def __repr__(self):
        return f"5' UTR: ({self.id}, {self.ensembl_id}) {self.start}-{self.end}"


    def __init__(self, start, end, strand, assembly_id, ensembl_id, refseq_id, 
                chess_id, velia_id, attrs={}):
       
        super(UtrFive, self).__init__(start, end, strand, assembly_id, 'utr_five')
        self.ensembl_id = ensembl_id
        self.refseq_id = refseq_id
        self.chess_id = chess_id
        self.velia_id = velia_id
        self.attrs = attrs


class Cds(SequenceRegion):
    """Represents a coding sequence (CDS) annotation.

    Stores information about protein-coding regions including various
    protein identifiers.

    Attributes:
        id (int): Primary key, inherited from SequenceRegion
        ensembl_id (str): Ensembl CDS ID
        refseq_id (str): RefSeq CDS ID
        chess_id (str): CHESS database ID
        openprot_id (str): OpenProt database ID
        velia_id (str): Velia database ID
        ccds_id (str): CCDS database ID
        ensembl_protein_id (str): Ensembl protein ID
        refseq_protein_id (str): RefSeq protein ID
        attrs (JSON): Additional attributes
    """
    """
    """
    __tablename__ = 'cds'

    id = Column(Integer, ForeignKey('sequence_region.id'), primary_key=True)
    ensembl_id  = Column(String(60))
    refseq_id = Column(String(60))
    chess_id = Column(String(60))
    openprot_id = Column(String(60))
    velia_id = Column(String)
    ccds_id  = Column(String(60))
    ensembl_protein_id = Column(String(30))
    refseq_protein_id = Column(String(30))
    attrs = Column(JSON)

    __mapper_args__ = { 'polymorphic_identity': 'cds' }

    def __repr__(self):
        return f"CDS:({self.id}) ({self.assembly_id}, {self.strand}) {self.start}-{self.end}"

    def __init__(self, start, end, strand, assembly_id, ensembl_id, refseq_id, 
                 chess_id, velia_id, ccds_id, ensembl_protein_id,
                 refseq_protein_id, openprot_id='', attrs={}):
        
        super(Cds, self).__init__(start, end, strand, assembly_id, 'cds')
        self.ensembl_id = ensembl_id
        self.refseq_id = refseq_id
        self.chess_id = chess_id
        self.velia_id = velia_id
        self.ccds_id = ccds_id
        self.ensembl_protein_id = ensembl_protein_id
        self.refseq_protein_id = refseq_protein_id
        self.openprot_id = openprot_id
        self.attrs = attrs


class Orf(Base):
    """Represents an Open Reading Frame (ORF).

    Stores information about potential protein-coding regions including
    sequence data and various identifiers.

    Attributes:
        id (int): Primary key
        start (int): Start position
        end (int): End position
        strand (str): Strand ('+' or '-')
        assembly_id (int): Foreign key to Assembly
        block_sizes (str): Sizes of ORF blocks
        chrom_starts (str): Chromosome start positions
        phases (str): Phase information
        reading_frames (str): Reading frame information
        orf_idx (str): Unique ORF index
        orf_idx_str (str): Human-readable ORF index
        secondary_orf_id (str): Secondary ORF identifier
        aa_seq (str): Amino acid sequence
        nt_seq (str): Nucleotide sequence
        ensembl_protein_id (str): Ensembl protein ID
        refseq_protein_id (str): RefSeq protein ID
        uniprot_id (str): UniProt ID
        attrs (JSON): Additional attributes
        vtx_id (str): VTX identifier
    """
    __tablename__ = 'orf'

    id = Column(Integer, primary_key=True, autoincrement=True)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    strand = Column(String(1), nullable=False)
    assembly_id = Column(Integer, ForeignKey('assembly.id'))
    block_sizes = Column(String)
    chrom_starts = Column(String)
    phases = Column(String)
    reading_frames = Column(String)

    # f'{orf_start}_{orf_end}_{row.strand}_{row.assembly_id}_{block_sizes}_{chrom_starts}_{exon_frames}'
    orf_idx = Column(String(130))
    orf_idx_str = Column(String)

    secondary_orf_id = Column(String)
    aa_seq = Column(String)
    nt_seq = Column(String)
    ensembl_protein_id = Column(String(30))
    refseq_protein_id = Column(String(30))
    uniprot_id = Column(String(30))
    attrs = Column(JSON)
    vtx_id = Column(String)

    assembly = relationship('Assembly', primaryjoin='Orf.assembly_id==Assembly.id')
    xref = relationship('OrfXref', primaryjoin='Orf.id==OrfXref.orf_id')
    
    __table_args__ = (UniqueConstraint('orf_idx'), {})

    def __repr__(self):
        return f"ORF: {self.id}, {self.vtx_id}, {self.secondary_orf_id},"\
               f"({self.strand}){self.assembly.ucsc_style_name}:{self.start}-{self.end}, {self.orf_idx_str}"

    def __init__(self, start, end, strand, assembly_id, block_sizes, chrom_starts, phases, reading_frames, orf_idx, orf_idx_str, secondary_orf_id, 
                 aa_seq, nt_seq, ensembl_protein_id, refseq_protein_id, uniprot_id, attrs={}, vtx_id='', id=None):
        
        if id:
            self.id = id
        self.start = start
        self.end = end
        self.strand = strand
        self.assembly_id = assembly_id
        self.block_sizes = block_sizes
        self.chrom_starts = chrom_starts
        self.phases = phases
        self.reading_frames = reading_frames
        self.orf_idx = orf_idx
        self.orf_idx_str = orf_idx_str

        self.secondary_orf_id = secondary_orf_id
        self.aa_seq = aa_seq
        self.nt_seq = nt_seq
        self.ensembl_protein_id = ensembl_protein_id
        self.refseq_protein_id = refseq_protein_id
        self.uniprot_id = uniprot_id
        self.attrs = attrs
        self.vtx_id = vtx_id


class Protein(Base):
    """Represents a protein sequence.

    Stores protein sequence data and various database identifiers.

    Attributes:
        id (int): Primary key
        aa_seq (str): Amino acid sequence
        aa_idx (str): Unique amino acid sequence index
        ensembl_protein_id (str): Ensembl protein ID
        refseq_protein_id (str): RefSeq protein ID
        openprot_id (str): OpenProt database ID
        uniprot_id (str): UniProt ID
        velia_id (str): Velia database ID
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'protein'

    id = Column(Integer, primary_key=True, autoincrement=True)
    aa_seq = Column(String)
    aa_idx = Column(String)
    ensembl_protein_id = Column(String(30))
    refseq_protein_id = Column(String(30))
    openprot_id = Column(String(30))
    uniprot_id = Column(String(30))
    velia_id = Column(String(30))
    attrs = Column(JSON)

    __table_args__ = (UniqueConstraint('aa_idx'), {})

    def __repr__(self):
        return f"Protein (#{self.id}, {self.uniprot_id}, {self.ensembl_protein_id})"

    def __init__(self, aa_seq, aa_idx, ensembl_protein_id, refseq_protein_id, uniprot_id, openprot_id='', velia_id=''):
        self.aa_seq = aa_seq
        self.aa_idx = aa_idx
        self.ensembl_protein_id = ensembl_protein_id
        self.refseq_protein_id = refseq_protein_id
        self.uniprot_id = uniprot_id
        self.openprot_id = openprot_id
        self.velia_id = velia_id


class ProteinOrf(Base):
    """Junction table linking Protein and ORF tables.

    Attributes:
        protein_id (int): Foreign key to Protein
        orf_id (int): Foreign key to Orf
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'protein_orf'

    protein_id = Column(Integer, ForeignKey('protein.id'), primary_key=True)
    orf_id = Column(Integer, ForeignKey('orf.id'), primary_key=True)
    attrs = Column(JSON)

    def __init__(self, protein_id, orf_id, attrs={}):
        self.protein_id = protein_id
        self.orf_id = orf_id
        self.attrs = attrs


class TranscriptUtrFive(Base):
    """Junction table linking Transcript and UtrFive tables.

    Attributes:
        id (int): Primary key
        transcript_id (int): Foreign key to Transcript
        utr_five_id (int): Foreign key to UtrFive
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'transcript_utr_five'

    id = Column(Integer, primary_key=True, autoincrement=True)
    transcript_id = Column(Integer, ForeignKey('transcript.id'))
    utr_five_id = Column(Integer, ForeignKey('utr_five.id'))
    attrs = Column(JSON)

    __table_args__ = (UniqueConstraint('transcript_id', 'utr_five_id'), {})

    def __init__(self, transcript_id, utr_five_id, attrs={}):
        self.transcript_id = transcript_id
        self.utr_five_id = utr_five_id
        self.attrs = attrs


class TranscriptUtrThree(Base):
    """Junction table linking Transcript and UtrThree tables.

    Attributes:
        id (int): Primary key
        transcript_id (int): Foreign key to Transcript
        utr_three_id (int): Foreign key to UtrThree
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'transcript_utr_three'

    id = Column(Integer, primary_key=True, autoincrement=True)
    transcript_id = Column(Integer, ForeignKey('transcript.id'))
    utr_three_id = Column(Integer, ForeignKey('utr_three.id'))
    attrs = Column(JSON)

    __table_args__ = (UniqueConstraint('transcript_id', 'utr_three_id'), {})

    def __init__(self, transcript_id, utr_three_id, attrs={}):
        self.transcript_id = transcript_id
        self.utr_three_id = utr_three_id
        self.attrs = attrs


class TranscriptExon(Base):
    """Junction table linking Transcript and Exon tables.

    Attributes:
        id (int): Primary key
        transcript_id (int): Foreign key to Transcript
        exon_id (int): Foreign key to Exon
        exon_number (int): Order of exon in transcript
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'transcript_exon'

    id = Column(Integer, primary_key=True, autoincrement=True)
    transcript_id = Column(Integer, ForeignKey('transcript.id'))
    exon_id = Column(Integer, ForeignKey('exon.id'))
    exon_number = Column(Integer)
    attrs = Column(JSON)

    __table_args__ = (UniqueConstraint('transcript_id', 'exon_id', 'exon_number'), {})

    def __init__(self, transcript_id, exon_id, exon_number, attrs={}):
        self.transcript_id = transcript_id
        self.exon_id = exon_id
        self.exon_number = exon_number
        self.attrs = attrs


class TranscriptOrf(Base):
    """Junction table linking Transcript and ORF tables.

    Attributes:
        transcript_id (int): Foreign key to Transcript
        orf_id (int): Foreign key to Orf
        evidence_tag (str): Evidence for the association
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'transcript_orf'

    transcript_id = Column(Integer, ForeignKey('transcript.id'), primary_key=True)
    orf_id = Column(Integer, ForeignKey('orf.id'), primary_key=True)
    evidence_tag = Column(String)
    attrs = Column(JSON)

    def __init__(self, transcript_id, orf_id, evidence_tag):
        self.transcript_id = transcript_id
        self.orf_id = orf_id
        self.evidence_tag = evidence_tag


class ExonCds(Base):
    """Junction table linking Exon and CDS tables.

    Attributes:
        exon_id (int): Foreign key to Exon
        cds_id (int): Foreign key to Cds
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'exon_cds'

    exon_id = Column(Integer, ForeignKey('exon.id'), primary_key=True)
    cds_id = Column(Integer, ForeignKey('cds.id'), primary_key=True)
    attrs = Column(JSON)

    def __init__(self, exon_id, cds_id):
        self.exon_id = exon_id
        self.cds_id = cds_id


class CdsOrf(Base):
    """Junction table linking CDS and ORF tables.

    Attributes:
        cds_id (int): Foreign key to Cds
        orf_id (int): Foreign key to Orf
        cds_number (int): Order of CDS in ORF
        phase (int): Phase of the CDS
        reading_frame (int): Reading frame
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'cds_orf'

    cds_id = Column(Integer, ForeignKey('cds.id'), primary_key=True)
    orf_id = Column(Integer, ForeignKey('orf.id'), primary_key=True)
    cds_number = Column(Integer)
    phase = Column(Integer)
    reading_frame = Column(Integer)
    attrs = Column(JSON)

    def __init__(self, cds_id, orf_id, cds_number, phase, reading_frame=0):
        self.cds_id = cds_id
        self.orf_id = orf_id
        self.cds_number = cds_number
        self.phase = phase
        self.reading_frame = reading_frame


class Dataset(Base):
    """Represents a dataset in the database.

    Stores information about different types of datasets.

    Attributes:
        id (int): Primary key
        name (str): Dataset name
        description (str): Dataset description
        type (str): Dataset type
        attrs (JSON): Additional attributes
    """
    __tablename__ = 'dataset'

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String)
    description = Column(String)
    type = Column(String(40))
    attrs = Column(JSON)

    __mapper_args__ = {'polymorphic_identity': 'dataset',
                       'polymorphic_on': type}

    __table_args__ = (UniqueConstraint('name',),{})

    def __repr__(self):
        return f"Data Set ({self.id}): {self.name}"

    def __init__(self, name, description, type, attrs={}):
        self.name = name
        self.description = description
        self.type = type
        self.attrs = attrs


class OrfXref(Base):
    """Cross-reference table for ORFs.

    Links ORFs to external database identifiers.

    Attributes:
        id (int): Primary key
        orf_id (int): Foreign key to Orf
        xref (str): External reference identifier
        type (str): Type of cross-reference
        attrs (JSON): Additional attributes
        orf_dataset_id (int): Foreign key to source Dataset
        xref_dataset_id (int): Foreign key to target Dataset
    """
    __tablename__ = 'orf_xref'

    id = Column(Integer, primary_key=True, autoincrement=True)
    orf_id = Column(Integer, ForeignKey('orf.id'))
    xref = Column(String)
    type = Column(String(30))
    attrs = Column(JSON)
    
    orf_dataset_id = Column(Integer, ForeignKey('dataset.id'))
    xref_dataset_id = Column(Integer, ForeignKey('dataset.id'))

    orf = relationship('Orf', primaryjoin = orf_id == Orf.id, viewonly=True)
    orf_data_source = relationship('Dataset', primaryjoin = orf_dataset_id == Dataset.id)
    xref_data_source = relationship('Dataset', primaryjoin = xref_dataset_id == Dataset.id)

    __table_args__ = (UniqueConstraint('orf_id', 'xref', 'xref_dataset_id'),{})

    def __repr__(self):
        return f"xref:({self.orf_id}, {str(self.orf_data_source)}) <-> ({self.xref}, {self.xref_data_source})"

    def __init__(self, orf_id, xref, type, orf_dataset_id, xref_dataset_id, attrs={}):
        self.orf_id = orf_id
        self.xref = xref
        self.type = type
        self.orf_dataset_id = orf_dataset_id
        self.xref_dataset_id = xref_dataset_id
        self.attrs = attrs


class TranscriptXref(Base):
    """Cross-reference table for Transcripts.

    Links Transcripts to external database identifiers.

    Attributes:
        id (int): Primary key
        transcript_id (int): Foreign key to Transcript
        xref (str): External reference identifier
        type (str): Type of cross-reference
        attrs (JSON): Additional attributes
        transcript_dataset_id (int): Foreign key to source Dataset
        xref_dataset_id (int): Foreign key to target Dataset
    """
    __tablename__ = 'transcript_xref'

    id = Column(Integer, primary_key=True, autoincrement=True)
    transcript_id = Column(Integer, ForeignKey('transcript.id'))
    xref = Column(String)
    type = Column(String(30))
    attrs = Column(JSON)
    
    transcript_dataset_id = Column(Integer, ForeignKey('dataset.id'))
    xref_dataset_id = Column(Integer, ForeignKey('dataset.id'))

    transcript = relationship('Transcript', primaryjoin = transcript_id == Transcript.id, viewonly=True)
    transcript_data_source = relationship('Dataset', primaryjoin = transcript_dataset_id == Dataset.id)
    xref_data_source = relationship('Dataset', primaryjoin = xref_dataset_id == Dataset.id)

    __table_args__ = (UniqueConstraint('transcript_id', 'xref', 'xref_dataset_id'),{})

    def __repr__(self):
        return f"xref:({self.transcript_id}, {str(self.transcript_data_source)}) <-> ({self.xref}, {self.xref_data_source})"

    def __init__(self, transcript_id, xref, type, transcript_dataset_id, xref_dataset_id, attrs={}):
        self.transcript_id = transcript_id
        self.xref = xref
        self.type = type
        self.transcript_dataset_id = transcript_dataset_id
        self.xref_dataset_id = xref_dataset_id
        self.attrs = attrs


class SequenceRegionXref(Base):
    """Cross-reference table for SequenceRegions.

    Links SequenceRegions to external database identifiers.

    Attributes:
        id (int): Primary key
        sequence_region_id (int): Foreign key to SequenceRegion
        xref (str): External reference identifier
        type (str): Type of cross-reference
        attrs (JSON): Additional attributes
        sequence_region_dataset_id (int): Foreign key to source Dataset
        xref_dataset_id (int): Foreign key to target Dataset
    """
    __tablename__ = 'sequence_region_xref'

    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence_region_id = Column(Integer, ForeignKey('sequence_region.id'))
    xref = Column(String)
    type = Column(String(30))
    attrs = Column(JSON)
    
    sequence_region_dataset_id = Column(Integer, ForeignKey('dataset.id'))
    xref_dataset_id = Column(Integer, ForeignKey('dataset.id'))

    sequence_region_data_source = relationship('Dataset', primaryjoin = sequence_region_dataset_id == Dataset.id)
    xref_data_source = relationship('Dataset', primaryjoin = xref_dataset_id == Dataset.id)

    __table_args__ = (UniqueConstraint('sequence_region_id', 'xref', 'xref_dataset_id'),{})

    def __repr__(self):
        return f"SequenceRegionID: {self.sequence_region_id} => {self.type}: {self.xref}, {self.xref_data_source}"

    def __init__(self, sequence_region_id, xref, type, sequence_region_dataset_id, xref_dataset_id, attrs={}):
        self.sequence_region_id = sequence_region_id
        self.xref = xref
        self.type = type
        self.sequence_region_dataset_id = sequence_region_dataset_id
        self.xref_dataset_id = xref_dataset_id
        self.attrs = attrs


class ProteinXref(Base):
    """Cross-reference table for Proteins.

    Links Proteins to external database identifiers.

    Attributes:
        id (int): Primary key
        protein_id (int): Foreign key to Protein
        xref (str): External reference identifier
        type (str): Type of cross-reference
        attrs (JSON): Additional attributes
        protein_dataset_id (int): Foreign key to source Dataset
        xref_dataset_id (int): Foreign key to target Dataset
    """
    __tablename__ = 'protein_xref'

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey('protein.id'))
    xref = Column(String)
    type = Column(String(30))
    attrs = Column(JSON)
    
    protein_dataset_id = Column(Integer, ForeignKey('dataset.id'))
    xref_dataset_id = Column(Integer, ForeignKey('dataset.id'))

    protein_data_source = relationship('Dataset', primaryjoin = protein_dataset_id == Dataset.id)
    xref_data_source = relationship('Dataset', primaryjoin = xref_dataset_id == Dataset.id)

    __table_args__ = (UniqueConstraint('protein_id', 'xref', 'xref_dataset_id'),{})

    def __repr__(self):
        return f"xref:({self.protein_id}, {str(self.protein_data_source)}) <-> ({self.xref}, {self.xref_data_source})"

    def __init__(self, protein_id, xref, type, protein_dataset_id, xref_dataset_id, attrs={}):
        self.protein_id = protein_id
        self.xref = xref
        self.type = type
        self.protein_dataset_id = protein_dataset_id
        self.xref_dataset_id = xref_dataset_id
        self.attrs = attrs


class _Session(_SA_Session):
    """Custom SQLAlchemy session class for the Velia database.

    Extends the standard SQLAlchemy Session with additional convenience methods:
    - get_or_create: Get existing or create new objects
    - upsert: Perform PostgreSQL upsert operations

    Notes:
        - Automatically sets search_path to settings.schema
    """

    def __init__(self, *args, **kwargs):
        super(_Session, self).__init__(*args, **kwargs)
        self.get_or_create = MethodType(get_or_create, self)
        self.upsert = MethodType(upsert, self)
        #self.search_by_xref = MethodType(search_by_xref, self)


    def __repr__(self):
        return "VeliaDB session %d" % (self.__hash__())


def get_or_create(session, class_type, **kwargs):
    """Gets or creates an object in the database based on unique constraints.

    Args:
        session: The database session object
        class_type: The SQLAlchemy model class to query/create
        **kwargs: Keyword arguments containing the column values to filter/create with

    Returns:
        The existing or newly created database object

    Notes:
        - Requires the model class to have a UniqueConstraint defined
        - Will check inherited class constraints if the model inherits from another
        - Creates a new object if no existing one is found matching the constraints
    """

    for constraint in list(class_type.__table_args__):
        if constraint.__class__.__name__ == 'UniqueConstraint':
            unique_cols = constraint.columns.keys()

    inherited_result = True
    if '__mapper_args__' in class_type.__dict__ and 'inherits' in class_type.__mapper_args__:
        inherited_class_type = class_type.__mapper_args__['inherits']
        for constraint in list(inherited_class_type.__table_args__):
            if constraint.__class__.__name__ == 'UniqueConstraint':
                inherited_unique_cols = constraint.columns.keys()

        try: inherited_result = session.query(inherited_class_type).filter_by(**{k: kwargs[k] for k in inherited_unique_cols}).first()
        except: None

    result = session.query(class_type).filter_by(**{k: kwargs[k] for k in unique_cols}).first()

    if not result or not inherited_result:
        result = class_type(**kwargs)
        session.add(result)
        session.commit()

    return result


def update(session, object, **kwargs):
    """Updates an existing database object with new values.

    Args:
        session: The database session object
        object: The database object to update
        **kwargs: Keyword arguments containing the column values to update

    Returns:
        The updated database object

    Notes:
        - Must pass in the actual object to update
        - Updates all attributes specified in kwargs
        - Commits changes to the database
    """
    #result = session.query(class_type).filter_by(**kwargs).first()
    #result = session.query(class_type).filter_by(name=kwargs['name']).first()
    #if result is None: return

    for key,value in kwargs.iteritems():
        setattr(object,key,value)
    session.add(object)
    session.commit()

    return object


def upsert(session, class_type, **kwargs):
    """Performs an upsert (insert or update) operation using PostgreSQL ON CONFLICT.

    Args:
        session: The database session object
        class_type: The SQLAlchemy model class to upsert
        **kwargs: Keyword arguments containing the column values

    Returns:
        The upserted database object

    Notes:
        - Uses the first unique constraint as the conflict target
        - Updates all columns on conflict
        - Executes as a single database operation
    """
    for constraint in list(class_type.__table_args__):
        if constraint.__class__.__name__ == 'UniqueConstraint':
            unique_cols = constraint.columns.keys()
            
    stmt = insert(class_type).values([kwargs,])
    
    stmt = stmt.on_conflict_do_update(
        index_elements=[eval(f'{class_type.__name__}.{unique_cols[0]}')], set_={k: eval(f'stmt.excluded.{k}', {}, {'stmt': stmt}) for k in kwargs.keys()}
    )
    
    session.execute(stmt)
    session.commit()
    
    result = session.query(class_type).filter_by(**{k: kwargs[k] for k in unique_cols}).first()
    
    return result


Session = sessionmaker(bind=engine, class_=_Session)

