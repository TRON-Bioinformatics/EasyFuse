"""
Tests for fusion validation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.src.breakpoint import Breakpoint
from bin.fusionannotation.src.cds import CDS
from bin.fusionannotation.src.fusion_transcript import FusionTranscript, get_involved_features
from bin.fusionannotation.src.transcript import Transcript

class TestFusionTranscript(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.wt1 = Transcript("ENST00000456", "", "GeneA", "", "")
        self.wt2 = Transcript("ENST00000457", "", "GeneB", "", "")
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("1", 300, "+")
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)


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
        bp1 = Breakpoint("1", 2000, "+")
        bp2 = Breakpoint("1", 2500, "+")
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.determine_fusion_type(1000), "cis_near")


    def test_get_fusion_type_cis_far(self):
        """Test case for type definition"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("1", 300, "+")
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.determine_fusion_type(50), "cis_far")


    def test_get_fusion_type_cis_trans(self):
        """Test case for type definition"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("1", 100, "+")
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.determine_fusion_type(50), "cis_trans")


    def test_get_fusion_type_cis_inv(self):
        """Test case for type definition"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("1", 300, "-")
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.determine_fusion_type(50), "cis_inv")


    def test_get_fusion_type_trans(self):
        """Test case for type definition"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("2", 300, "+")
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.determine_fusion_type(100), "trans")


    def test_get_fusion_type_trans_inv(self):
        """Test case for type definition"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "-")
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        self.assertEqual(self.fusion_transcript.determine_fusion_type(100), "trans_inv")


    def test_get_combined_boundary_trans_inv(self):
        """Test case for boundary calculation"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "-")
        bp1.exon_boundary = "right_boundary"
        bp2.exon_boundary = "right_boundary"
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_boundary()
        self.assertEqual(result, "both")


    def test_get_combined_boundary_trans(self):
        """Test case for boundary calculation"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        bp1.exon_boundary = "right_boundary"
        bp2.exon_boundary = "right_boundary"
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_boundary()
        self.assertEqual(result, "no_match")


    def get_combined_frame_no_frame(self):
        """Test case for frame calculation"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        self.wt1.frame = "-1"
        self.wt2.frame = "-1"
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_frame()
        self.assertEqual(result, "no_frame")


    def get_combined_frame_neo_frame(self):
        """Test case for frame calculation"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        self.wt1.frame = "2"
        self.wt2.frame = "-1"
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_frame()
        self.assertEqual(result, "neo_frame")


    def get_combined_frame_in_frame(self):
        """Test case for frame calculation"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        self.wt1.frame = "2"
        self.wt2.frame = "2"
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_frame()
        self.assertEqual(result, "in_frame")


    def get_combined_frame_out_frame(self):
        """Test case for frame calculation"""
        bp1 = Breakpoint("1", 200, "+")
        bp2 = Breakpoint("3", 300, "+")
        self.wt1.frame = "1"
        self.wt2.frame = "2"
        self.fusion_transcript = FusionTranscript(self.wt1, self.wt2, bp1, bp2)
        result = self.fusion_transcript.get_combined_frame()
        self.assertEqual(result, "out_frame")


    def get_involved_features_cds_first_partner(self):
        """Test case for exon selection"""
        cds = [
            CDS(cds_id='cds_1', start=30157, stop=30179, frame=0, transcript_id='transcript_1')
        ]
        bp = Breakpoint("21", 30157, "-")
        actual_result = get_involved_features(cds, bp, False)
        expected_result = [
            CDS(cds_id='cds_1', start=30157, stop=30179, frame=0, transcript_id='transcript_1')
        ]
        self.assertEqual(actual_result, expected_result)


    def get_involved_features_cds_second_partner(self):
        """Test case for exon selection"""
        cds = [
            CDS(cds_id='cds_2', start=248000, stop=248030, frame=0, transcript_id='transcript_2')
        ]
        bp = Breakpoint("7", 248000, "+")
        actual_result = get_involved_features(cds, bp, True)
        expected_result = [
            CDS(cds_id='cds_2', start=248000, stop=248030, frame=0, transcript_id='transcript_2')
        ]
        self.assertEqual(actual_result, expected_result)
