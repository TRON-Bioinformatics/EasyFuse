"""
Tests for fusion validation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.src.breakpoint import Breakpoint
from bin.fusionannotation.src.fusion_transcript import FusionTranscript, get_involved_exons
from bin.fusionannotation.src.transcript import Transcript

class TestFusionTranscript(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("1", 300, "+")
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)


    def test_get_bpid(self):
        """Test case for breakpoint ID generation"""
        self.assertEqual(self.fusion_transcript.get_bpid(), "1:200:+_1:300:+")


    def test_get_ftid(self):
        """Test case for fusion transcript ID generation"""
        ftid = "GeneA_1:200:+_ENST00000456_GeneB_1:300:+_ENST00000457"
        self.assertEqual(self.fusion_transcript.get_ftid(), ftid)


    def test_get_fusion_gene_name(self):
        """Test case for fusion gene name generation"""
        self.assertEqual(self.fusion_transcript.get_fusion_gene_name(), "GeneA_GeneB")


    def test_get_fusion_type_cis_near(self):
        """Test case for type definition"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 2000, "+")
        bp2 = Breakpoint("1", 2500, "+")
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.get_fusion_type(1000), "cis_near")


    def test_get_fusion_type_cis_far(self):
        """Test case for type definition"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("1", 300, "+")
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.get_fusion_type(50), "cis_far")


    def test_get_fusion_type_cis_trans(self):
        """Test case for type definition"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("1", 100, "+")
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.get_fusion_type(50), "cis_trans")


    def test_get_fusion_type_cis_inv(self):
        """Test case for type definition"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("1", 300, "-")
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.get_fusion_type(50), "cis_inv")


    def test_get_fusion_type_trans(self):
        """Test case for type definition"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("2", 300, "+")
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.get_fusion_type(100), "trans")


    def test_get_fusion_type_trans_inv(self):
        """Test case for type definition"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "-")
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.get_fusion_type(100), "trans_inv")


    def test_get_combined_boundary_trans_inv(self):
        """Test case for boundary calculation"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "-")
        bp1.exon_boundary = "right_boundary"
        bp2.exon_boundary = "right_boundary"
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_boundary()
        self.assertEqual(result, "both")


    def test_get_combined_boundary_trans(self):
        """Test case for boundary calculation"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        bp1.exon_boundary = "right_boundary"
        bp2.exon_boundary = "right_boundary"
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_boundary()
        self.assertEqual(result, "no_match")


    def get_combined_frame_no_frame(self):
        """Test case for frame calculation"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        wt1.frame = "-1"
        wt2.frame = "-1"
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_frame()
        self.assertEqual(result, "no_frame")


    def get_combined_frame_neo_frame(self):
        """Test case for frame calculation"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        wt1.frame = "2"
        wt2.frame = "-1"
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_frame()
        self.assertEqual(result, "neo_frame")


    def get_combined_frame_in_frame(self):
        """Test case for frame calculation"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        wt1.frame = "2"
        wt2.frame = "2"
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_frame()
        self.assertEqual(result, "in_frame")


    def get_combined_frame_out_frame(self):
        """Test case for frame calculation"""
        wt1 = Transcript("transcript:ENST00000456", "", "GeneA", "", "")
        wt2 = Transcript("transcript:ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        wt1.frame = "1"
        wt2.frame = "2"
        self.fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_frame()
        self.assertEqual(result, "out_frame")
