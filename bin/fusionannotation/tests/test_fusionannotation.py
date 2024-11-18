"""
Tests for fusion annotation module.
"""

import unittest

# pylint: disable=E0401
import bin.fusionannotation.fusionannotation # type: ignore


class TestFusionAnnotation(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.input_fusions = "tests/fusionannotation/Detected_Fusions.csv"

    def test_run(self):
        """Test case for BP loading."""
        pass
