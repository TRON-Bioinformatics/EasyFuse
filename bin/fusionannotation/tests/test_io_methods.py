"""
Tests for IO methods module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.src.breakpoint import Breakpoint
from bin.fusionannotation.src.io_methods import load_detected_fusions # type: ignore

class TestIOMethods(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.input_fusions = "tests/fusionannotation/input_data/Detected_Fusions.csv"

    def test_load_detected_fusions(self):
        """Test case for BP loading."""
        res = {
            '21:41494380:-_7:13935843:-': (
                Breakpoint('21', 41494380, '-'), Breakpoint('7', 13935843, '-')
            ),
            '22:29287134:+_12:50814280:+': (
                Breakpoint('22', 29287134, '+'), Breakpoint('12', 50814280, '+')
            ),
            '7:92003235:+_7:140787584:-': (
                Breakpoint('7', 92003235, '+'), Breakpoint('7', 140787584, '-')
            ),
            '12:11869969:+_15:87940753:-': (
                Breakpoint('12', 11869969, '+'), Breakpoint('15', 87940753, '-')
            ),
            '8:42968214:+_10:43116584:+': (
                Breakpoint('8', 42968214, '+'), Breakpoint('10', 43116584, '+')
            ),
            '22:29287134:+_11:128807180:+': (
                Breakpoint('22', 29287134, '+'), Breakpoint('11', 128807180, '+')
            ),
            '19:15254152:-_15:34347969:+': (
                Breakpoint('19', 15254152, '-'), Breakpoint('15', 34347969, '+')
            ),
            '5:150404680:-_6:117324415:-': (
                Breakpoint('5', 150404680, '-'), Breakpoint('6', 117324415, '-')
            ),
            '21:41494375:-_7:13935838:-': (
                Breakpoint('21', 41494375, '-'), Breakpoint('7', 13935838, '-')
            ),
            '2:42301394:+_2:29223584:-': (
                Breakpoint('2', 42301394, '+'), Breakpoint('2', 29223584, '-')
            )
        }
        self.assertEqual(load_detected_fusions(self.input_fusions), res)
