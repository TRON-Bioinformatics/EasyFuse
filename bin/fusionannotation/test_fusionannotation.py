"""
Tests for fusion annotation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.io_methods import breakpoints_to_dict

class TestFusionAnnotation(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.input_fusions = "tests/fusionannotation/Detected_Fusions.csv"

    def test_breakpoints_to_dict(self):
        """Test case for BP loading."""
        res = {
            '21:41494380:-_7:13935843:-': ('21:41494380:-', '7:13935843:-'), 
            '22:29287134:+_12:50814280:+': ('22:29287134:+', '12:50814280:+'), 
            '7:92003235:+_7:140787584:-': ('7:92003235:+', '7:140787584:-'), 
            '12:11869969:+_15:87940753:-': ('12:11869969:+', '15:87940753:-'), 
            '8:42968214:+_10:43116584:+': ('8:42968214:+', '10:43116584:+'), 
            '22:29287134:+_11:128807180:+': ('22:29287134:+', '11:128807180:+'), 
            '19:15254152:-_15:34347969:+': ('19:15254152:-', '15:34347969:+'), 
            '5:150404680:-_6:117324415:-': ('5:150404680:-', '6:117324415:-'), 
            '21:41494375:-_7:13935838:-': ('21:41494375:-', '7:13935838:-'), 
            '2:42301394:+_2:29223584:-': ('2:42301394:+', '2:29223584:-')
        }
        self.assertEqual(breakpoints_to_dict(self.input_fusions), res)
