"""
Tests for fusion validation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.fusion_validation import define_type
from bin.fusionannotation.fusion_validation import get_fusion_feature_list

class TestIOMethods(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.input_fusions = "tests/fusionannotation/Detected_Fusions.csv"

    def test_define_type(self):
        """Test case for type definition"""
        pass


    def test_get_fusion_feature_list(self):
        """Test case for fusion feature list generation"""
        pass
