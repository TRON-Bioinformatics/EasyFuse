"""
Tests for exon validation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.exon_validation import check_exon_overlap
from bin.fusionannotation.exon_validation import get_wt_codings
from bin.fusionannotation.exon_validation import get_bp_overlapping_features
from bin.fusionannotation.exon_validation import get_parents


class TestExonValidation(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.input_fusions = "tests/fusionannotation/Detected_Fusions.csv"

    def test_check_exon_overlap(self):
        """Test case for checking exon overlap"""
        pass


    def test_get_wt_codings(self):
        """Test case for wt codings"""
        pass


    def test_get_bp_overlapping_features(self):
        """Test case for overlapping features"""
        pass


    def test_get_parents(self):
        """Test case for parents determination"""
        pass