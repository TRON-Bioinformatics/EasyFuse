"""
Tests for fusion validation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.src.breakpoint import Breakpoint
from bin.fusionannotation.src.cds import CDS
from bin.fusionannotation.src.exon import Exon

class TestBreakpoint(unittest.TestCase):
    """
    Provides unit tests for breakpoint module.
    """
    def setUp(self):
        self.bp = Breakpoint("1", 100, "+")


    def test_get_boundary_left(self):
        """Test case for left boundary check."""
        exon = Exon("1", 100, 200, "ENST00000456")
        self.assertEqual(self.bp.get_boundary(exon), "left_boundary")


    def test_get_boundary_right(self):
        """Test case for right boundary check."""
        exon = Exon("2", 1, 100, "ENST00000456")
        self.assertEqual(self.bp.get_boundary(exon), "right_boundary")


    def test_get_boundary_within(self):
        """Test case for boundary lying within the exon."""
        exon = Exon("3", 50, 150, "ENST00000456")
        self.assertEqual(self.bp.get_boundary(exon), "within")


    def test_get_boundary_outside(self):
        """Test case for boundary lying outside the exon."""
        exon = Exon("4", 101, 200, "ENST00000456")
        self.assertEqual(self.bp.get_boundary(exon), "outside")


    def test_get_boundary_none(self):
        """Test case for non-existent exon."""
        exon = None
        self.assertEqual(self.bp.get_boundary(exon), "NA")


    def test_get_frame(self):
        """Test case for frame calculation."""
        cds = [
            CDS("CDS:ENSP00000504602_11", 1, 2, 1,"transcript:ENST00000679263"),
            CDS("CDS:ENSP00000504238_9", 4, 5, 2, "transcript:ENST00000679181"),
        ]
        result = self.bp.get_frame(cds)
        frame = (-1, -1)
        self.assertEqual(result, frame)


    def test_get_frame_none(self):
        """Test case for non-existent CDS."""
        cds = []
        result = self.bp.get_frame(cds)
        frame = (-1, -1)
        self.assertEqual(result, frame)


    def test_get_frame_within(self):
        """Test case for frame calculation within CDS."""
        cds = [
            CDS("CDS:ENSP00000504602_11", 50, 100, 1,"transcript:ENST00000679263"),
            CDS("CDS:ENSP00000504238_9", 75, 125, 1, "transcript:ENST00000679181"),
            CDS("CDS:ENSP00000504705_11", 3, 4, 2, "transcript:ENST00000676973"),
        ]
        result = self.bp.get_frame(cds)
        frame = (-1, -1)
        self.assertEqual(result, frame)
