"""
Tests for exon validation module.
"""

import unittest

# pylint: disable=E0401
from src.fusionannotation.cds import CDS
from src.fusionannotation.exon import Exon
from src.fusionannotation.feature_validation import get_exon_cds_overlap
from src.fusionannotation.feature_validation import filter_cds_by_exons


class TestFeatureValidation(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """

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
        result, _ = get_exon_cds_overlap(exons, cds)
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


if __name__ == '__main__':
    unittest.main()
