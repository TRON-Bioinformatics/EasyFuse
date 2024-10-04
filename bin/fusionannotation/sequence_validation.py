"""
This module provides methods for sequence operations.
"""

def get_stranded_seq(sequence, strand):
    """
    Return the reverse complement of a sequence 
    if the strand is negative and the unchanged sequence otherwise.
    """
    if strand == "-":
        return sequence.reverse_complement()
    return sequence
