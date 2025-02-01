"""
This module provides utility functions for handling and manipulating data structures such as lists and dictionaries.
It also includes a function for writing JSONL results.
"""
import csv
import gzip
import json
import logging
from typing import Iterable, Dict, List, Callable, Any, Type, Union, Collection, Optional, IO

logger = logging.getLogger(__name__)


def count(items: Iterable) -> Dict:
    """
    Counts the occurrence of each item in the iterable.

    Args:
        items (Iterable): The iterable to count items from.

    Returns:
        Dict: A dictionary where the keys are the items and the values are the counts.
    """
    f = {}
    for item in items:
        if item not in f:
            f[item] = 1
        else:
            f[item] += 1
    return f


def sort_dict(unsorted_dict: Dict) -> Dict:
    """
    Sorts a dictionary by its keys.

    Args:
        unsorted_dict (Dict): The dictionary to sort.

    Returns:
        Dict: The sorted dictionary.
    """
    return {k: unsorted_dict[k] for k in sorted(unsorted_dict.keys())}


def get_iterator_chunk(iterator: Iterable[object], chunk_size: int) -> List:
    """
    Gets a chunk of a specified size from an iterator.

    Args:
        iterator (Iterable[object]): The iterator to get a chunk from.
        chunk_size (int): The size of the chunk to get.

    Returns:
        List: The chunk from the iterator.
    """
    chunk = []

    while len(chunk) < chunk_size:
        try:
            # noinspection PyTypeChecker
            chunk.append(next(iterator))
        except StopIteration:
            break

    return chunk


def generate_jsonl_result_writer(outfile: Any) -> Callable:
    """
    Generates a function that writes JSONL results to an output file.

    Args:
        outfile (Any): The output file to write to.

    Returns:
        Callable: The function that writes JSONL results.
    """

    def write_orfs(orfset: Collection):
        logger.info(f'Writing {len(orfset)} ORFs to {outfile.filename} ...')
        for this_orf in orfset:
            outfile.write(json.dumps(this_orf.to_dict()) + '\n')

    return write_orfs


def get_dict_chunk(this_dict: Dict, keys: Iterable) -> Dict:
    """
    Gets a chunk of a dictionary using a set of keys.

    Args:
        this_dict (Dict): The dictionary to get a chunk from.
        keys (Iterable): The keys to use for getting the chunk.

    Returns:
        Dict: The chunk from the dictionary.
    """
    return {k: this_dict[k] for k in keys}


def freq(items_to_count: Iterable) -> dict:
    counts = count(items_to_count)
    total_count = sum(counts.values())
    return {item: this_count / total_count for item, this_count in counts.items()}


def extract_numeric_vector(field_string: str, seps: List[str] = None,
                           data_type: Type[Union[int, float]
                                           ] = lambda x: int(float(x)),
                           remove_dangling_separators: bool = True) -> List[Union[int, float]]:
    """
    Extracts a numeric vector from a string.
    Args:
        field_string: The string to extract from.
        seps: The separators to use. Defaults to None.
        data_type: The data type to convert to. Defaults to int.
        remove_dangling_separators: Whether to remove dangling separators. Defaults to False.
    Returns:
        The extracted numeric vector.
    """
    if seps is None:
        seps = [';', ',']
    for sep in seps:
        try:
            if remove_dangling_separators:
                parsed_vals = [
                    data_type(element) for element in field_string.strip(sep).split(sep)]
            else:
                parsed_vals = [data_type(element)
                               for element in field_string.split(sep)]

        except ValueError:
            continue
        else:
            return parsed_vals
    raise ValueError(
        f'Could not find a valid way to parse {field_string} using separators {seps} with data type {data_type}.')


def encode_numeric_vector(vector_elements: List[Union[int, float]], sep: str = ';') -> str:
    """
    Encodes a numeric vector into a string.
    Args:
        vector_elements: The elements of the vector.
        sep: The separator to use. Defaults to ';'.
    Returns:
        The encoded string.
    """
    return sep.join([str(e) for e in vector_elements])


class LazyCsvWriter:
    """
    A class for lazily writing rows to a CSV file.
    Field names can be provided at initialization or will be inferred from the first row of data.

    Attributes:
        output_fpath (str): The file path to the output CSV file.
        fieldnames (List[str]): The list of fieldnames for the CSV file. Defaults to None.
        writer (csv.DictWriter): The CSV writer object. Initialized as None.
    """

    def __init__(self, output_fpath: str, name: Optional[str] = '', fieldnames: Optional[List[str]] = None) -> None:
        """
        Initializes the LazyCsvWriter object.

        Args:
            output_fpath (str): The file path to the output CSV file.
            fieldnames (Optional[List[str]]): The list of fieldnames for the CSV file. Defaults to None.
        """
        self.name = name
        self.output_fpath: str = output_fpath
        self.fieldnames: Optional[List[str]] = fieldnames
        self.file = None
        self.file: Optional[IO] = None
        self.writer: Optional[csv.DictWriter] = None

    def initialize_writer(self) -> None:
        """
        Initializes the CSV writer object and writes the header to the CSV file.
        """
        if self.output_fpath.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        self.file = opener(self.output_fpath, 'wt')
        self.writer = csv.DictWriter(self.file, fieldnames=self.fieldnames)
        self.writer.writeheader()
        logger.info('Writing %s to %s with fields %s ...',
                    self.name, self.output_fpath, ','.join(self.fieldnames))

    def __del__(self) -> None:
        if self.file:
            self.file.close()

    def writerows(self, rows: List[Dict[Any, Any]]) -> None:
        """
        Writes a list of rows to the CSV file.

        Args:
            rows (List[Dict[Any, Any]]): The rows to write to the CSV file.
        """
        if not len(rows):
            return
        if not self.fieldnames:
            self.fieldnames = sorted(rows[0].keys())
        if not self.writer:
            self.initialize_writer()
        self.writer.writerows(rows)
