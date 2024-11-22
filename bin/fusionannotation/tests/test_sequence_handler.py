"""
Tests for IO methods module.
"""

import unittest

# pylint: disable=E0401
from Bio.Seq import Seq # type: ignore

from src.fusionannotation.breakpoint import Breakpoint
from src.fusionannotation.sequence_handler import calc_hash
from src.fusionannotation.sequence_handler import concatenate_seqs
from src.fusionannotation.sequence_handler import get_context_sequence
from src.fusionannotation.sequence_handler import get_fusion_transcript_sequence
from src.fusionannotation.sequence_handler import get_peptide_sequence
from src.fusionannotation.sequence_handler import get_stranded_seq

class TestSequenceHandler(unittest.TestCase):
    """
    Provides unit tests for the sequence handler module.
    """


    def test_calc_hash(self):
        """Test case for hash calculation"""
        seq = Seq('ACGTTTA')
        result = calc_hash(seq)
        self.assertEqual(result, 'a40d472865950d7c')


    def test_get_stranded_seq_positive_strand(self):
        """Test case for strandedness validation"""
        self.assertEqual(get_stranded_seq(Seq("ATCG"), "+"), Seq("ATCG"))


    def test_get_stranded_seq_negative_strand(self):
        """Test case for strandedness validation"""
        self.assertEqual(get_stranded_seq(Seq("ATCG"), "-"), Seq("CGAT"))


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
        """Test case for peptide sequence validation"""
        cds_seq = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
        bp = Breakpoint("21", 30157, "+")
        result = get_peptide_sequence(cds_seq, bp, 0)
        self.assertEqual(result, Seq("MAIVMGR*KGAR*"))


    def test_get_fusion_transcript_sequence_in_frame(self):
        """Test case for fusion transcript sequence validation"""
        ft1_cds_seq = Seq('ACGT')
        ft2_cds_seq = Seq('ATTA')
        ft2_exon_seq = Seq('GTCG')

        result = get_fusion_transcript_sequence(
            ft1_cds_seq,
            ft2_cds_seq,
            ft2_exon_seq,
            "-",
            "+",
            "in_frame"
        )
        self.assertEqual(result, Seq("ACGTATTA"))


    def test_get_fusion_transcript_sequence_out_frame(self):
        """Test case for fusion transcript sequence validation"""
        ft1_cds_seq = Seq('ACGT')
        ft2_cds_seq = Seq('ATTA')
        ft2_exon_seq = Seq('GTCG')

        result = get_fusion_transcript_sequence(
            ft1_cds_seq,
            ft2_cds_seq,
            ft2_exon_seq,
            "-",
            "+",
            "neo_frame"
        )
        self.assertEqual(result, Seq("ACGTGTCG"))


if __name__ == '__main__':
    unittest.main()
