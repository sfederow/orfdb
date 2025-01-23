"""Module to implement ORM for the Velia Database"""

from typing import Dict, Any, Optional
from orfdb.settings import db_connection_string

from sqlalchemy.orm import (
    sessionmaker, 
    relationship, 
    DeclarativeBase,
    Mapped,
    mapped_column
)
from sqlalchemy import (
    create_engine,
    ForeignKey, 
    String, 
    Integer, 
    JSON,
    UniqueConstraint, 
    select
)
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.pool import NullPool
from types import MethodType


def create_db_engine(connection_string: str):
    """Create SQLAlchemy engine with 2.0 style configuration.
    
    Args:
        connection_string: Database connection URL
    """
    engine = create_engine(
        connection_string,
        future=True,
        poolclass=NullPool,
    )
    return engine


def get_session_maker(engine):
    """Create a session factory configured for SQLAlchemy 2.0.
    
    Args:
        engine: SQLAlchemy engine instance
    """
    return sessionmaker(
        bind=engine,
        expire_on_commit=False,
        future=True
    )


engine = create_db_engine(db_connection_string)
Session = get_session_maker(engine)


class Base(DeclarativeBase):
    """Base class for all SQLAlchemy models"""
    pass


class Assembly(Base):
    """Represents a genomic assembly in the database."""
    
    __tablename__ = 'assembly'

    id: Mapped[int] = mapped_column(primary_key=True)
    genbank_accession: Mapped[Optional[str]] = mapped_column(String(200))
    refseq_accession: Mapped[Optional[str]] = mapped_column(String(100))
    ucsc_style_name: Mapped[Optional[str]] = mapped_column(String(50))
    sequence_role: Mapped[Optional[str]] = mapped_column(String(200))
    assembly_unit: Mapped[Optional[str]] = mapped_column(String(200))
    assigned_molecule: Mapped[Optional[str]] = mapped_column(String(200))
    sequence_length: Mapped[int] = mapped_column(Integer)
    genome_accession: Mapped[Optional[str]] = mapped_column(String(200))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    __table_args__ = (
        UniqueConstraint('genbank_accession', 'refseq_accession'),
    )

    def __repr__(self) -> str:
        return (f"Assembly (#{self.id}): {self.genome_accession}, "
                f"{self.ucsc_style_name}, {self.genbank_accession}, "
                f"{self.sequence_length}")


class SequenceRegion(Base):
    """Base class for genomic sequence regions."""
    
    __tablename__ = 'sequence_region'

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    assembly_id: Mapped[int] = mapped_column(ForeignKey('assembly.id'))
    start: Mapped[int] = mapped_column(Integer)
    end: Mapped[int] = mapped_column(Integer)
    strand: Mapped[str] = mapped_column(String(1))
    type: Mapped[str] = mapped_column(String(30))

    __table_args__ = (
        UniqueConstraint('start', 'end', 'strand', 'assembly_id', 'type'),
    )

    __mapper_args__ = {
        'polymorphic_identity': 'sequence_region',
        'polymorphic_on': type
    }

    def __repr__(self) -> str:
        return f"SequenceRegion: {self.start}-{self.end} ({self.strand})"

    def __repr__dict__(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "start": self.start,
            "end": self.end,
            "strand": self.strand
        }


class Gene(SequenceRegion):
    """Represents a gene annotation."""
    
    __tablename__ = 'gene'

    id: Mapped[int] = mapped_column(ForeignKey('sequence_region.id'), primary_key=True)
    hgnc_id: Mapped[Optional[str]] = mapped_column(String(30))
    hgnc_name: Mapped[Optional[str]] = mapped_column(String(30))
    ensembl_id: Mapped[Optional[str]] = mapped_column(String(30))
    refseq_id: Mapped[Optional[str]] = mapped_column(String(30))
    chess_id: Mapped[Optional[str]] = mapped_column(String(30))
    velia_id: Mapped[Optional[str]] = mapped_column(String(30))
    gene_type: Mapped[Optional[str]] = mapped_column(String(200))
    long_name: Mapped[Optional[str]] = mapped_column(String(100))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    transcripts: Mapped[list["Transcript"]] = relationship(
        'Transcript', 
        primaryjoin='Gene.id==Transcript.gene_id'
    )

    __mapper_args__ = {'polymorphic_identity': 'gene'}

    def __repr__(self) -> str:
        return (f"Gene: ({self.hgnc_id}, {self.ensembl_id}, {self.hgnc_name})"
                f"{self.start}-{self.end} ({self.strand})")


class Transcript(Base):
    """Represents a transcript annotation."""
    
    __tablename__ = 'transcript'

    id: Mapped[int] = mapped_column(primary_key=True)
    start: Mapped[int] = mapped_column(Integer)
    end: Mapped[int] = mapped_column(Integer)
    strand: Mapped[str] = mapped_column(String(1))
    assembly_id: Mapped[int] = mapped_column(ForeignKey('assembly.id'))
    block_sizes: Mapped[Optional[str]] = mapped_column(String)
    chrom_starts: Mapped[Optional[str]] = mapped_column(String)
    transcript_idx: Mapped[str] = mapped_column(String(130))
    transcript_idx_str: Mapped[Optional[str]] = mapped_column(String)

    gene_id: Mapped[int] = mapped_column(ForeignKey('gene.id'))
    ensembl_id: Mapped[Optional[str]] = mapped_column(String(30))
    refseq_id: Mapped[Optional[str]] = mapped_column(String(30))
    chess_id: Mapped[Optional[str]] = mapped_column(String(30))
    velia_id: Mapped[Optional[str]] = mapped_column(String(30))
    support_level: Mapped[Optional[str]] = mapped_column(String)
    transcript_type: Mapped[Optional[str]] = mapped_column(String)
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    gene: Mapped["Gene"] = relationship('Gene', viewonly=True)
    exons: Mapped[list["Exon"]] = relationship('Exon', secondary='transcript_exon', viewonly=True)

    __table_args__ = (UniqueConstraint('transcript_idx'),)

    def __repr__(self) -> str:
        return (f"Transcript: ({self.id}, {self.ensembl_id}, {self.refseq_id})"
                f"{self.strand}:{self.start}-{self.end}")


class Exon(SequenceRegion):
    """Represents an exon annotation."""
    
    __tablename__ = 'exon'

    id: Mapped[int] = mapped_column(ForeignKey('sequence_region.id'), primary_key=True)
    ensembl_id: Mapped[Optional[str]] = mapped_column(String(30))
    refseq_id: Mapped[Optional[str]] = mapped_column(String(30))
    chess_id: Mapped[Optional[str]] = mapped_column(String(30))
    velia_id: Mapped[Optional[str]] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    transcripts: Mapped[list["Transcript"]] = relationship(
        'Transcript', 
        secondary='transcript_exon', 
        viewonly=True
    )

    __mapper_args__ = {'polymorphic_identity': 'exon'}

    def __repr__(self) -> str:
        return f"Exon: ({self.id}, {self.ensembl_id}) {self.start}-{self.end}"


class UtrThree(SequenceRegion):
    """Represents a 3' UTR annotation."""
    
    __tablename__ = 'utr_three'

    id: Mapped[int] = mapped_column(ForeignKey('sequence_region.id'), primary_key=True)
    ensembl_id: Mapped[Optional[str]] = mapped_column(String(30))
    refseq_id: Mapped[Optional[str]] = mapped_column(String(30))
    chess_id: Mapped[Optional[str]] = mapped_column(String(30))
    velia_id: Mapped[Optional[str]] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    __mapper_args__ = {'polymorphic_identity': 'utr_three'}

    def __repr__(self) -> str:
        return f"3' UTR: ({self.id}, {self.ensembl_id}) {self.start}-{self.end}"


class UtrFive(SequenceRegion):
    """Represents a 5' UTR annotation."""
    
    __tablename__ = 'utr_five'

    id: Mapped[int] = mapped_column(ForeignKey('sequence_region.id'), primary_key=True)
    ensembl_id: Mapped[Optional[str]] = mapped_column(String(30))
    refseq_id: Mapped[Optional[str]] = mapped_column(String(30))
    chess_id: Mapped[Optional[str]] = mapped_column(String(30))
    velia_id: Mapped[Optional[str]] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    __mapper_args__ = {'polymorphic_identity': 'utr_five'}

    def __repr__(self) -> str:
        return f"5' UTR: ({self.id}, {self.ensembl_id}) {self.start}-{self.end}"


class Cds(SequenceRegion):
    """Represents a coding sequence (CDS) annotation."""
    
    __tablename__ = 'cds'

    id: Mapped[int] = mapped_column(ForeignKey('sequence_region.id'), primary_key=True)
    ensembl_id: Mapped[Optional[str]] = mapped_column(String(60))
    refseq_id: Mapped[Optional[str]] = mapped_column(String(60))
    chess_id: Mapped[Optional[str]] = mapped_column(String(60))
    openprot_id: Mapped[Optional[str]] = mapped_column(String(60))
    velia_id: Mapped[Optional[str]] = mapped_column(String)
    ccds_id: Mapped[Optional[str]] = mapped_column(String(60))
    ensembl_protein_id: Mapped[Optional[str]] = mapped_column(String(30))
    refseq_protein_id: Mapped[Optional[str]] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    __mapper_args__ = {'polymorphic_identity': 'cds'}

    def __repr__(self) -> str:
        return f"CDS:({self.id}) ({self.assembly_id}, {self.strand}) {self.start}-{self.end}"


class Orf(Base):
    """Represents an Open Reading Frame (ORF)."""
    
    __tablename__ = 'orf'

    id: Mapped[int] = mapped_column(primary_key=True)
    start: Mapped[int] = mapped_column(Integer)
    end: Mapped[int] = mapped_column(Integer)
    strand: Mapped[str] = mapped_column(String(1))
    assembly_id: Mapped[int] = mapped_column(ForeignKey('assembly.id'))
    block_sizes: Mapped[Optional[str]] = mapped_column(String)
    chrom_starts: Mapped[Optional[str]] = mapped_column(String)
    phases: Mapped[Optional[str]] = mapped_column(String)
    reading_frames: Mapped[Optional[str]] = mapped_column(String)
    orf_idx: Mapped[str] = mapped_column(String(130))
    orf_idx_str: Mapped[Optional[str]] = mapped_column(String)
    secondary_orf_id: Mapped[Optional[str]] = mapped_column(String)
    aa_seq: Mapped[Optional[str]] = mapped_column(String)
    nt_seq: Mapped[Optional[str]] = mapped_column(String)
    ensembl_protein_id: Mapped[Optional[str]] = mapped_column(String(30))
    refseq_protein_id: Mapped[Optional[str]] = mapped_column(String(30))
    uniprot_id: Mapped[Optional[str]] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)
    vtx_id: Mapped[Optional[str]] = mapped_column(String)

    assembly: Mapped["Assembly"] = relationship('Assembly')
    xref: Mapped[list["OrfXref"]] = relationship('OrfXref')

    __table_args__ = (UniqueConstraint('orf_idx'),)

    def __repr__(self) -> str:
        return (f"ORF: {self.id}, {self.vtx_id}, {self.secondary_orf_id},"
                f"({self.strand}){self.assembly.ucsc_style_name}:{self.start}-{self.end}, "
                f"{self.orf_idx_str}")


class Protein(Base):
    """Represents a protein sequence."""
    
    __tablename__ = 'protein'

    id: Mapped[int] = mapped_column(primary_key=True)
    aa_seq: Mapped[str] = mapped_column(String)
    aa_idx: Mapped[str] = mapped_column(String)
    ensembl_protein_id: Mapped[Optional[str]] = mapped_column(String(30))
    refseq_protein_id: Mapped[Optional[str]] = mapped_column(String(30))
    openprot_id: Mapped[Optional[str]] = mapped_column(String(30))
    uniprot_id: Mapped[Optional[str]] = mapped_column(String(30))
    velia_id: Mapped[Optional[str]] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    __table_args__ = (UniqueConstraint('aa_idx'),)

    def __repr__(self) -> str:
        return f"Protein (#{self.id}, {self.uniprot_id}, {self.ensembl_protein_id})"


class ProteinOrf(Base):
    """Junction table linking Protein and ORF tables."""
    
    __tablename__ = 'protein_orf'

    protein_id: Mapped[int] = mapped_column(ForeignKey('protein.id'), primary_key=True)
    orf_id: Mapped[int] = mapped_column(ForeignKey('orf.id'), primary_key=True)
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)


class TranscriptUtrFive(Base):
    """Junction table linking Transcript and UtrFive tables."""
    
    __tablename__ = 'transcript_utr_five'

    id: Mapped[int] = mapped_column(primary_key=True)
    transcript_id: Mapped[int] = mapped_column(ForeignKey('transcript.id'))
    utr_five_id: Mapped[int] = mapped_column(ForeignKey('utr_five.id'))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    __table_args__ = (UniqueConstraint('transcript_id', 'utr_five_id'),)


class TranscriptUtrThree(Base):
    """Junction table linking Transcript and UtrThree tables."""
    
    __tablename__ = 'transcript_utr_three'

    id: Mapped[int] = mapped_column(primary_key=True)
    transcript_id: Mapped[int] = mapped_column(ForeignKey('transcript.id'))
    utr_three_id: Mapped[int] = mapped_column(ForeignKey('utr_three.id'))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    __table_args__ = (UniqueConstraint('transcript_id', 'utr_three_id'),)


class TranscriptExon(Base):
    """Junction table linking Transcript and Exon tables."""
    
    __tablename__ = 'transcript_exon'

    id: Mapped[int] = mapped_column(primary_key=True)
    transcript_id: Mapped[int] = mapped_column(ForeignKey('transcript.id'))
    exon_id: Mapped[int] = mapped_column(ForeignKey('exon.id'))
    exon_number: Mapped[int]
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    __table_args__ = (UniqueConstraint('transcript_id', 'exon_id', 'exon_number'),)


class TranscriptOrf(Base):
    """Junction table linking Transcript and ORF tables."""
    
    __tablename__ = 'transcript_orf'

    transcript_id: Mapped[int] = mapped_column(ForeignKey('transcript.id'), primary_key=True)
    orf_id: Mapped[int] = mapped_column(ForeignKey('orf.id'), primary_key=True)
    evidence_tag: Mapped[Optional[str]] = mapped_column(String)
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)


class ExonCds(Base):
    """Junction table linking Exon and CDS tables."""
    
    __tablename__ = 'exon_cds'

    exon_id: Mapped[int] = mapped_column(ForeignKey('exon.id'), primary_key=True)
    cds_id: Mapped[int] = mapped_column(ForeignKey('cds.id'), primary_key=True)
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)


class CdsOrf(Base):
    """Junction table linking CDS and ORF tables."""
    
    __tablename__ = 'cds_orf'

    cds_id: Mapped[int] = mapped_column(ForeignKey('cds.id'), primary_key=True)
    orf_id: Mapped[int] = mapped_column(ForeignKey('orf.id'), primary_key=True)
    cds_number: Mapped[int] = mapped_column(Integer)
    phase: Mapped[int] = mapped_column(Integer)
    reading_frame: Mapped[int] = mapped_column(Integer)
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)


class Dataset(Base):
    """Represents a dataset in the database."""
    
    __tablename__ = 'dataset'

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String)
    description: Mapped[str] = mapped_column(String)
    type: Mapped[str] = mapped_column(String(40))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)

    __mapper_args__ = {
        'polymorphic_identity': 'dataset',
        'polymorphic_on': type
    }

    __table_args__ = (UniqueConstraint('name'),)

    def __repr__(self) -> str:
        return f"Data Set ({self.id}): {self.name}"


class OrfXref(Base):
    """Cross-reference table for ORFs."""
    
    __tablename__ = 'orf_xref'

    id: Mapped[int] = mapped_column(primary_key=True)
    orf_id: Mapped[int] = mapped_column(ForeignKey('orf.id'))
    xref: Mapped[str] = mapped_column(String)
    type: Mapped[str] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)
    
    orf_dataset_id: Mapped[int] = mapped_column(ForeignKey('dataset.id'))
    xref_dataset_id: Mapped[int] = mapped_column(ForeignKey('dataset.id'))

    orf: Mapped["Orf"] = relationship('Orf', viewonly=True)
    orf_data_source: Mapped["Dataset"] = relationship(
        'Dataset', 
        primaryjoin='OrfXref.orf_dataset_id==Dataset.id'
    )
    xref_data_source: Mapped["Dataset"] = relationship(
        'Dataset', 
        primaryjoin='OrfXref.xref_dataset_id==Dataset.id'
    )

    __table_args__ = (UniqueConstraint('orf_id', 'xref', 'xref_dataset_id'),)

    def __repr__(self) -> str:
        return (f"xref:({self.orf_id}, {str(self.orf_data_source)}) <-> "
                f"({self.xref}, {self.xref_data_source})")
    

class TranscriptXref(Base):
    """Cross-reference table for Transcripts."""
    
    __tablename__ = 'transcript_xref'

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    transcript_id: Mapped[int] = mapped_column(ForeignKey('transcript.id'))
    xref: Mapped[str] = mapped_column(String)
    type: Mapped[str] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)
    
    transcript_dataset_id: Mapped[int] = mapped_column(ForeignKey('dataset.id'))
    xref_dataset_id: Mapped[int] = mapped_column(ForeignKey('dataset.id'))

    transcript: Mapped["Transcript"] = relationship(
        'Transcript', 
        primaryjoin='TranscriptXref.transcript_id==Transcript.id', 
        viewonly=True
    )
    transcript_data_source: Mapped["Dataset"] = relationship(
        'Dataset', 
        primaryjoin='TranscriptXref.transcript_dataset_id==Dataset.id'
    )
    xref_data_source: Mapped["Dataset"] = relationship(
        'Dataset', 
        primaryjoin='TranscriptXref.xref_dataset_id==Dataset.id'
    )

    __table_args__ = (UniqueConstraint('transcript_id', 'xref', 'xref_dataset_id'),)

    def __repr__(self) -> str:
        return f"xref:({self.transcript_id}, {str(self.transcript_data_source)}) <-> ({self.xref}, {self.xref_data_source})"


class SequenceRegionXref(Base):
    """Cross-reference table for SequenceRegions."""
    
    __tablename__ = 'sequence_region_xref'

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    sequence_region_id: Mapped[int] = mapped_column(ForeignKey('sequence_region.id'))
    xref: Mapped[str] = mapped_column(String)
    type: Mapped[str] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)
    
    sequence_region_dataset_id: Mapped[int] = mapped_column(ForeignKey('dataset.id'))
    xref_dataset_id: Mapped[int] = mapped_column(ForeignKey('dataset.id'))

    sequence_region_data_source: Mapped["Dataset"] = relationship(
        'Dataset', 
        primaryjoin='SequenceRegionXref.sequence_region_dataset_id==Dataset.id'
    )
    xref_data_source: Mapped["Dataset"] = relationship(
        'Dataset', 
        primaryjoin='SequenceRegionXref.xref_dataset_id==Dataset.id'
    )

    __table_args__ = (UniqueConstraint('sequence_region_id', 'xref', 'xref_dataset_id'),)

    def __repr__(self) -> str:
        return f"SequenceRegionID: {self.sequence_region_id} => {self.type}: {self.xref}, {self.xref_data_source}"


class ProteinXref(Base):
    """Cross-reference table for Proteins."""
    
    __tablename__ = 'protein_xref'

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    protein_id: Mapped[int] = mapped_column(ForeignKey('protein.id'))
    xref: Mapped[str] = mapped_column(String)
    type: Mapped[str] = mapped_column(String(30))
    attrs: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSON)
    
    protein_dataset_id: Mapped[int] = mapped_column(ForeignKey('dataset.id'))
    xref_dataset_id: Mapped[int] = mapped_column(ForeignKey('dataset.id'))

    protein_data_source: Mapped["Dataset"] = relationship(
        'Dataset', 
        primaryjoin='ProteinXref.protein_dataset_id==Dataset.id'
    )
    xref_data_source: Mapped["Dataset"] = relationship(
        'Dataset', 
        primaryjoin='ProteinXref.xref_dataset_id==Dataset.id'
    )

    __table_args__ = (UniqueConstraint('protein_id', 'xref', 'xref_dataset_id'),)

    def __repr__(self) -> str:
        return f"xref:({self.protein_id}, {str(self.protein_data_source)}) <-> ({self.xref}, {self.xref_data_source})"


def upsert(session, class_type: type[Base], **kwargs: Any) -> Base:
    """Upsert (insert or update) a record using unique constraints.
    
    Args:
        session: SQLAlchemy session
        class_type: SQLAlchemy model class
        **kwargs: Column values to insert/update
        
    Returns:
        The inserted or updated model instance
        
    Note:
        Uses the first unique constraint found in the model's __table_args__
        for conflict resolution.
    """
    # Get unique constraint columns
    unique_cols = []
    for constraint in class_type.__table_args__:
        if isinstance(constraint, UniqueConstraint):
            unique_cols = constraint.columns.keys()
            break
            
    # Build the upsert statement
    stmt = insert(class_type).values([kwargs])
    
    # Create the on conflict do update clause
    stmt = stmt.on_conflict_do_update(
        index_elements=[getattr(class_type, col) for col in unique_cols],
        set_=kwargs
    )
    
    # Execute the statement
    session.execute(stmt)
    session.commit()
    
    # Query for the resulting record
    result = session.scalars(
        select(class_type).filter_by(
            **{k: kwargs[k] for k in unique_cols}
        )
    ).first()
    
    return result
