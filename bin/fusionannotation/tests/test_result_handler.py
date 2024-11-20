"""
Tests for result handler module.
"""

import unittest

# pylint: disable=E0401
from Bio.Seq import Seq # type: ignore

from bin.fusionannotation.src.exon import Exon
from bin.fusionannotation.src.fusion_transcript import FusionTranscript
from bin.fusionannotation.src.transcript import Transcript
from bin.fusionannotation.src.breakpoint import Breakpoint
from bin.fusionannotation.src.result_handler import generate_temp_dict
from bin.fusionannotation.src.result_handler import get_context_sequence
from bin.fusionannotation.src.result_handler import concatenate_seqs

class TestResultHandler(unittest.TestCase):
    """
    Provides unit tests for the sequence validation module.
    """

    def test_generate_temp_dict(self):
        """Test case for strandedness validation"""
        cds_seqs = {
            '21': [
                {
                    Exon(exon_id='exon_1', start=30157, stop=30179, transcript_id='transcript_1')
                },
                [
                    Seq('TTGAGGTGCACACTGTCCCGGAT')
                ]
            ],
            '7': [
                {
                    Exon(exon_id='exon_2', start=248000, stop=248030, transcript_id='transcript_2')
                },
                [
                    Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT')
                ]
            ]
        }
        exon = Exon(exon_id='exon_1', start=30157, stop=30179, transcript_id='transcript_1')
        tmp_dict = {
            exon: Seq('TTGAGGTGCACACTGTCCCGGAT')
        }
        result = generate_temp_dict(cds_seqs, '21')
        self.assertEqual(result, tmp_dict)


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
