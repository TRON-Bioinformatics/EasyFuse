"""
Tests for IO methods module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.src.sequence_validation import get_stranded_seq

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
