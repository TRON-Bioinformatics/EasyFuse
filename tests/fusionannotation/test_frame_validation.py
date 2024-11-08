"""
Tests for frame validation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.frame_validation import get_frame
from bin.fusionannotation.frame_validation import get_frame2


class TestFrameValidation(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.input_fusions = "tests/fusionannotation/Detected_Fusions.csv"

    def test_get_frame(self):
        """Test case for frame validation"""
        pass


    def test_get_frame2(self):
        """Test case for second frame validation"""
        pass
