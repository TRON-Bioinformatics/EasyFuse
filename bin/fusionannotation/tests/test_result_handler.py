"""
Tests for result handler module.
"""

import unittest

# pylint: disable=E0401
from Bio.Seq import Seq # type: ignore

from bin.fusionannotation.src.exon import Exon
from bin.fusionannotation.src.result_handler import generate_temp_dict


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
