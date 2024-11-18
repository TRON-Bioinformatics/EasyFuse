"""
Tests for exon validation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.src.cds import CDS
from bin.fusionannotation.src.exon import Exon
from bin.fusionannotation.src.feature_validation import check_exon_overlap
from bin.fusionannotation.src.feature_validation import get_exon_cds_overlap
from bin.fusionannotation.src.feature_validation import filter_cds_by_exons


class TestFeatureValidation(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """

    def test_check_exon_overlap_true(self):
        """Test case for checking exon overlap"""
        wt1_exons = [
            Exon("1", 1, 2, "transcript:ENST00000679263"),
            Exon("2", 4, 5, "transcript:ENST00000679181"),
            Exon("3", 7, 8, "transcript:ENST00000676973")
        ]
        wt2_exons = [
            Exon("1", 1, 2, "transcript:ENST00000679263"),
            Exon("2", 4, 5, "transcript:ENST00000679181"),
            Exon("3", 7, 8, "transcript:ENST00000676973")
        ]
        result = check_exon_overlap(wt1_exons, wt2_exons, "cis")
        self.assertTrue(result)


    def test_check_exon_overlap_false(self):
        """Test case for checking exon overlap"""
        wt1_exons = [
            Exon("1", 1, 2, "transcript:ENST00000679263"),
            Exon("2", 4, 5, "transcript:ENST00000679181"),
            Exon("3", 7, 8, "transcript:ENST00000676973")
        ]
        wt2_exons = [
            Exon("1", 9, 10, "transcript:ENST00000679263"),
            Exon("2", 11, 12, "transcript:ENST00000679181"),
            Exon("3", 13, 14, "transcript:ENST00000676973")
        ]
        result = check_exon_overlap(wt1_exons, wt2_exons, "cis")
        self.assertFalse(result)


    def test_check_exon_overlap_trans(self):
        """Test case for checking exon overlap"""
        wt1_exons = [
            Exon("1", 3, 4, "transcript:ENST00000679263")
        ]
        wt2_exons = [
            Exon("1", 1, 12, "transcript:ENST00000679263")
        ]
        result = check_exon_overlap(wt1_exons, wt2_exons, "trans")
        self.assertFalse(result)


    def test_get_exon_cds_overlap(self):
        """Test case for wt codings"""
        exons = [
            Exon("1", 1, 2, "transcript:ENST00000679263"),
            Exon("2", 4, 5, "transcript:ENST00000679181"),
            Exon("3", 7, 8, "transcript:ENST00000676973")
        ]
        cds = [
            CDS("1", 1, 3, 0, "transcript:ENST000006792"),
            CDS("2", 4, 5, 0, "transcript:ENST000006791"),
            CDS("3", 7, 8, 0, "transcript:ENST000006769")
        ]
        result, flagged_transcripts = get_exon_cds_overlap(exons, cds)
        expected = [
            None,
            CDS('2', 4, 5, 0, 'transcript:ENST000006791'),
            CDS('3', 7, 8, 0, 'transcript:ENST000006769')
        ]
        self.assertEqual(result, expected)


    def test_filter_cds_by_exons(self):
        """Test case for overlapping features"""
        exons = [
            Exon("1", 1, 2, "transcript:ENST00000679263"),
            Exon("2", 4, 5, "transcript:ENST00000679181"),
            Exon("3", 7, 8, "transcript:ENST00000676973")
        ]

        cds = [
            CDS("CDS:ENSP00000504602_11", 1, 2, 0, "transcript:ENST00000679263"),
            CDS("CDS:ENSP00000504238_9", 4, 5, 0, "transcript:ENST00000679181")
        ]
        result = filter_cds_by_exons(exons, cds)
        filtered_cds = [
            CDS('CDS:ENSP00000504602_11', 1, 2, 0, 'transcript:ENST00000679263'),
            CDS('CDS:ENSP00000504238_9', 4, 5, 0, 'transcript:ENST00000679181'),
            None
        ]
        self.assertEqual(result, filtered_cds)
