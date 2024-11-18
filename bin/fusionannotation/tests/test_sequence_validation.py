"""
Tests for IO methods module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.sequence_validation import get_stranded_seq

class TestSequenceValidation(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.input_fusions = "tests/fusionannotation/Detected_Fusions.csv"

    def test_get_stranded_seq(self):
        """Test case for strandedness validation"""
        pass
