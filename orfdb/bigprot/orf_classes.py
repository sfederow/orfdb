import gzip
import json
import os
import logging
from orfdb.bigprot import validation
import Bio.Seq
from typing import Dict, List, Union

logger = logging.getLogger(__name__)


class OrfBase:
    """
    Class representing an Open Reading Frame (ORF).
    The coordinates are pythonic, i.e. 0-based and half-open intervals (the end position is not in the sequence).
    To convert to GFF/Ensembl coordinates, add 1 to the start position and 2 to the end position.
    """
    FIELDS = [
        'codon_length',
        'nuc_length',
        'nuc_sequence',
        'parent_sequence_id',
        'phase',
        'start_codon',
        'start_pos',
        'end_pos',
        'stop_codon',
        'strand'
    ]

    def __init__(self, parent_sequence_id: str, start_pos: int, strand: str, codons: List[str]):
        """
        Initialize an OrfBase instance.
        :param parent_sequence_id: ID of the parent sequence.
        :param start_pos: Start position of the ORF.
        :param strand: Strand of the ORF.
        :param codons: List of codons in the ORF.
        """
        assert strand in {'+', '-'}, f'Invalid strand character {strand}!'
        self._num_to_codon: Dict[int, str] = {}
        self._num_to_codon = {}
        self._codon_to_num = {}
        self._curr_codon_num = 0

        self.parent_sequence_id = parent_sequence_id
        self.start_pos = start_pos
        validation.validate_strand(strand)
        self.strand = strand

        self._codons = []
        for codon in codons:
            self.add_codon(codon)

    @property
    def codons(self):
        """
        Get the codons in the ORF.
        :return: List of codons in the ORF.
        """
        return [self._num_to_codon[num] for num in self._codons]

    def add_codon(self, codon: str):
        """
        Add a codon to the ORF.
        :param codon: Codon to add.
        """
        if codon not in self._codon_to_num:
            self._curr_codon_num += 1
            self._num_to_codon[self._curr_codon_num] = codon
            self._codon_to_num[codon] = self._curr_codon_num

        self._codons.append(self._codon_to_num[codon])

    @property
    def start_codon(self):
        """
        Get the start codon of the ORF.
        :return: Start codon of the ORF, or None if the ORF has no codons.
        """
        if len(self._codons):
            return self._num_to_codon[self._codons[0]]
        else:
            return None

    @property
    def stop_codon(self):
        """
        Get the stop codon of the ORF.
        :return: Stop codon of the ORF, or None if the ORF has no codons.
        """
        if len(self._codons):
            return self._num_to_codon[self._codons[-1]]
        else:
            return None

    @property
    def codon_length(self):
        """
        Get the length of the ORF in codons.
        :return: Length of the ORF in codons.
        """
        return len(self._codons)

    @property
    def nuc_length(self):
        """
        Get the length of the ORF in nucleotides.
        :return: Length of the ORF in nucleotides.
        """
        return self.codon_length * 3

    @property
    def end_pos(self):
        """
        Get the end position of the ORF.
        :return: End position of the ORF.
        """
        return self.start_pos + self.nuc_length

    @property
    def nuc_sequence(self):
        """
        Get the nucleotide sequence of the ORF.
        :return: Nucleotide sequence of the ORF.
        """
        return ''.join(self.codons)

    @property
    def aa_sequence(self):
        """
        Get the amino acid sequence of the ORF.
        :return: Amino acid sequence of the ORF.
        """
        return Bio.Seq.translate(self.nuc_sequence, table='Standard')

    def __repr__(self):
        """
        Get a string representation of the ORF.
        :return: String representation of the ORF.
        """
        return f'ORF {self.parent_sequence_id} {self.strand} {self.start_pos}-{self.end_pos} {self.start_codon}-{self.stop_codon}'

    @property
    def unique_id(self):
        """
        Get a unique ID for the ORF.
        :return: Unique ID for the ORF.
        """
        return f'{self.parent_sequence_id}-{self.strand}-{self.start_pos}-{self.end_pos}'

    @property
    def description(self):
        """
        Get a description of the ORF.
        :return: Description of the ORF.
        """
        return f'ORF from {self.parent_sequence_id} {self.start_pos}-{self.end_pos} on {self.strand} with start codon {self.start_codon} and stop codon {self.stop_codon}'

    def to_dict(self):
        """
        Convert the ORF to a dictionary.
        :return: Dictionary representation of the ORF.
        """
        self_dict = {field: self.__getattribute__(
            field) for field in self.FIELDS}

        return self_dict

    @classmethod
    def from_dict(cls, attribute_dict: Dict[str, Union[str, int, List[str]]]):
        """
        Create an OrfBase instance from a dictionary.
        :param attribute_dict: Dictionary of attributes.
        :return: New OrfBase instance.
        """
        new_orf = cls(parent_sequence_id=attribute_dict['parent_sequence_id'],
                      start_pos=attribute_dict['start_pos'],
                      strand=attribute_dict['strand'],
                      codons=[attribute_dict['nuc_sequence'][pos:pos + 3] for pos in
                              range(0, len(attribute_dict['nuc_sequence']), 3)])
        # QC Checks
        assert new_orf.codon_length == attribute_dict['codon_length']

        return new_orf


class OrfSet:
    """
    Class representing a set of ORFs.
    """

    def __init__(self, orfs: List[OrfBase]):
        """
        Initialize an OrfSet instance.
        :param orfs: List of ORFs.
        """
        self.orfs = orfs

    def __len__(self):
        """
        Get the number of ORFs in the set.
        :return: Number of ORFs in the set.
        """
        return len(self.orfs)

    def to_dicts(self):
        """
        Convert the ORFs in the set to dictionaries.
        :return: List of dictionaries representing the ORFs in the set.
        """
        return [orf.to_dict() for orf in self.orfs]

    @classmethod
    def from_jsons(cls, jsons: List[Dict[str, Union[str, int, List[str]]]]):
        """
        Create an OrfSet instance from a list of JSON objects.
        :param jsons: List of JSON objects representing ORFs.
        :return: New OrfSet instance.
        """
        return cls([OrfBase.from_dict(json_dict) for json_dict in jsons])

    @classmethod
    def from_jsonl_file(cls, jsonl_fpath: str):
        """
        Create an OrfSet instance from a JSONL file.
        :param jsonl_fpath: Path to the JSONL file.
        :return: New OrfSet instance.
        """
        logger.info(f'Loading ORFs from {jsonl_fpath} ...')

        if os.path.splitext(jsonl_fpath)[-1] == '.gz':
            file_opener = gzip.open
        else:
            file_opener = open

        with file_opener(jsonl_fpath, 'rt') as jsonl_file:
            new_orfset = cls.from_jsons(
                [json.loads(line) for line in jsonl_file])

        logger.info(f'Loaded {len(new_orfset)} ORFs ...')
        return new_orfset
