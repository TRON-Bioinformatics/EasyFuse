"""
Tests for boundary validation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.boundary_validation import get_boundary
from bin.fusionannotation.boundary_validation import get_boundary2

class TestBoundaryValidation(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.input_fusions = "tests/fusionannotation/Detected_Fusions.csv"

    def test_get_boundary(self):
        """Test case for boundary determination"""
        pass


    def test_get_boundary2(self):
        """Test case for second boundary determination"""
        pass
