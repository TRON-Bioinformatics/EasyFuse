"""
Tests for IO methods module.
"""

import unittest

# pylint: disable=E0401
from Bio.Seq import Seq # type: ignore

from bin.fusionannotation.src.breakpoint import Breakpoint
from bin.fusionannotation.src.sequence_handler import calc_hash
from bin.fusionannotation.src.sequence_handler import concatenate_seqs
from bin.fusionannotation.src.sequence_handler import get_context_sequence
from bin.fusionannotation.src.sequence_handler import get_peptide_sequence
from bin.fusionannotation.src.sequence_handler import get_stranded_seq

class TestSequenceValidation(unittest.TestCase):
    """
    Provides unit tests for the sequence validation module.
    """


    def test_get_stranded_seq_positive_strand(self):
        """Test case for strandedness validation"""
        self.assertEqual(get_stranded_seq("ATCG", "+"), "ATCG")


    def test_get_stranded_seq_negative_strand(self):
        """Test case for strandedness validation"""
        self.assertEqual(get_stranded_seq("ATCG", "-"), "CGAT")


    def test_get_feature_transcripts(self):
        """Test case for strandedness validation"""
        seqs = [
            Seq('TTGAGGTGCACACTGTCCCGGAT'),
            Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT')
        ]
        concat_seq = Seq('TTGAGGTGCACACTGTCCCGGATTCTATCAAAATAAAAAATAGTGACAGCAAGT')
        result = concatenate_seqs(seqs)
        self.assertEqual(result, concat_seq)


    def test_get_feature_transcripts_empty(self):
        """Test case for strandedness validation"""
        seqs = []
        result = concatenate_seqs(seqs)
        self.assertEqual(result, Seq(''))


    def test_get_context_sequence(self):
        """Test case for context sequence validation"""
        ft1_seq = Seq('ACGTTTA')
        ft2_seq = Seq('GGATC')
        context_seq = Seq('TAAACGTGGATC')
        bp1 = Breakpoint("21", 30157, "-")
        bp2 = Breakpoint("7", 257882, "+")
        result = get_context_sequence(
            ft1_seq,
            ft2_seq,
            bp1.strand,
            bp2.strand
        )
        self.assertEqual(result, context_seq)


    def test_get_peptide_sequence(self):
        """Test case for context sequence validation"""
        cds_seq = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
        bp = Breakpoint("21", 30157, "+")
        result = get_peptide_sequence(cds_seq, bp, 0)
        self.assertEqual(result, Seq("MAIVMGR*KGAR*"))


    def test_calc_hash(self):
        """Test case for hash calculation"""
        seq = Seq('ACGTTTA')
        result = calc_hash(seq)
        self.assertEqual(result, 'a40d472865950d7c')
