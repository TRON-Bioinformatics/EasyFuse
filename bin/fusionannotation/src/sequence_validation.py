"""
This module provides methods for sequence operations.
"""

# pylint: disable=import-error
from Bio.Seq import Seq # type: ignore

def get_stranded_seq(sequence: Seq, strand: str) -> Seq:
    """
    Return the reverse complement of a sequence 
    if the strand is negative and the unchanged sequence otherwise.
    """
    if strand == "-":
        return sequence.reverse_complement()
    return sequence
