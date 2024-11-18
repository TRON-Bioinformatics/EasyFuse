"""
Tests for fusion validation module.
"""

import unittest

# pylint: disable=E0401
from bin.fusionannotation.src.breakpoint import Breakpoint # type: ignore
from bin.fusionannotation.src.cds import CDS
from bin.fusionannotation.src.exon import Exon
from bin.fusionannotation.src.gff_db_controller import GffDbController
from bin.fusionannotation.src.transcript import Transcript


class TestGffDbController(unittest.TestCase):
    """
    Provides unit tests for fusion annotation module.
    """
    def setUp(self):
        self.db = GffDbController("tests/fusionannotation/Homo_sapiens.GRCh38.110.gff3.db")


    def test_get_features_overlapping_position_exon(self):
        """Test case for extracting features overlapping the breakpoint"""
        # 21:41494380:-_7:13935843:- -> TMPRSS2_ETV1
        # TODO: create minimal test DB
        bp = Breakpoint("21", 41494380, "-")
        result = self.db.get_features_overlapping_position(bp, "exon")

        exons = [
            Exon("exon_1012058", 41494356, 41494578, "transcript:ENST00000679263"),
            Exon("exon_1012071", 41494356, 41494578, "transcript:ENST00000679181"),
            Exon("exon_1012084", 41494356, 41494578, "transcript:ENST00000676973"),
            Exon("exon_1012098", 41494356, 41494578, "transcript:ENST00000679054"),
            Exon("exon_1012112", 41494356, 41494578, "transcript:ENST00000678348"),
            Exon("exon_1012126", 41494356, 41494578, "transcript:ENST00000332149"),
            Exon("exon_1012138", 41494356, 41494578, "transcript:ENST00000679016"),
            Exon("exon_1012165", 41494356, 41494578, "transcript:ENST00000677680"),
            Exon("exon_1012179", 41494356, 41494578, "transcript:ENST00000678171"),
            Exon("exon_1012207", 41494356, 41494578, "transcript:ENST00000398585"),
            Exon("exon_1012221", 41494356, 41494578, "transcript:ENST00000458356"),
            Exon("exon_1012235", 41494356, 41494578, "transcript:ENST00000678617"),
            Exon("exon_1012247", 41494356, 41494578, "transcript:ENST00000424093"),
            Exon("exon_1012261", 41494356, 41494578, "transcript:ENST00000454499"),
            Exon("exon_1012271", 41494356, 41494578, "transcript:ENST00000455813")
        ]
        self.assertEqual(result, exons)


    def test_get_features_overlapping_position_cds(self):
        """Test case for extracting features overlapping the breakpoint"""
        # 21:41494380:-_7:13935843:- -> TMPRSS2_ETV1
        bp = Breakpoint("21", 41494380, "-")
        result = self.db.get_features_overlapping_position(bp, "CDS")

        cds = [
            CDS("CDS:ENSP00000504602_11", 41494356, 41494578, 0,"transcript:ENST00000679263"),
            CDS("CDS:ENSP00000504238_9", 41494356, 41494578, 0, "transcript:ENST00000679181"),
            CDS("CDS:ENSP00000504705_11", 41494356, 41494578, 0, "transcript:ENST00000676973"),
            CDS("CDS:ENSP00000502928_11", 41494356, 41494578, 0, "transcript:ENST00000679054"),
            CDS("CDS:ENSP00000503556_11", 41494356, 41494578, 0, "transcript:ENST00000678348"),
            CDS("CDS:ENSP00000330330_11", 41494356, 41494578, 0, "transcript:ENST00000332149"),
            CDS("CDS:ENSP00000504610_9", 41494356, 41494578, 0, "transcript:ENST00000679016"),
            CDS("CDS:ENSP00000504526_5", 41494356, 41494578, 0, "transcript:ENST00000677680"),
            CDS("CDS:ENSP00000503877_11", 41494356, 41494578, 0, "transcript:ENST00000678171"),
            CDS("CDS:ENSP00000381588_11", 41494356, 41494578, 0, "transcript:ENST00000398585"),
            CDS("CDS:ENSP00000391216_11", 41494356, 41494578, 0, "transcript:ENST00000458356"),
            CDS("CDS:ENSP00000397846_10", 41494356, 41494578, 0, "transcript:ENST00000424093"),
            CDS("CDS:ENSP00000389006_11", 41494356, 41494578, 0, "transcript:ENST00000454499"),
            CDS("CDS:ENSP00000391784", 41494356, 41494578, 0, "transcript:ENST00000455813")
        ]
        self.assertEqual(result, cds)


    def test_get_exons_from_transcript(self):
        """Test case for extracting exons from transcript"""
        transcript_id = "transcript:ENST00000679263"
        result = self.db.get_features_from_transcript(transcript_id, "exon")
        exons = [
            Exon('exon_1012047', 41464300, 41466153, ''),
            Exon('exon_1012048', 41467734, 41467886, ''),
            Exon('exon_1012049', 41468396, 41468538, ''),
            Exon('exon_1012050', 41470648, 41470743, ''),
            Exon('exon_1012051', 41471806, 41471981, ''),
            Exon('exon_1012052', 41473325, 41473496, ''),
            Exon('exon_1012053', 41476577, 41476620, ''),
            Exon('exon_1012054', 41479172, 41479282, ''),
            Exon('exon_1012055', 41480476, 41480602, ''),
            Exon('exon_1012056', 41488394, 41488513, ''),
            Exon('exon_1012057', 41489507, 41489593, ''),
            Exon('exon_1012058', 41494356, 41494578, ''),
            Exon('exon_1012059', 41498119, 41498189, ''),
            Exon('exon_1012060', 41511863, 41512108, '')
        ]
        self.assertEqual(result, exons)


    def test_get_cds_from_transcripts(self):
        """Test case for extracting cds from transcript"""
        transcript_id = "transcript:ENST00000679263"
        result = self.db.get_features_from_transcript(transcript_id, "CDS")

        cds = [
            CDS('CDS:ENSP00000504602', 41466142, 41466153, 0, ''),
            CDS('CDS:ENSP00000504602_1', 41467734, 41467886, 0, ''),
            CDS('CDS:ENSP00000504602_2', 41468396, 41468538, 2, ''),
            CDS('CDS:ENSP00000504602_3', 41470648, 41470743, 2, ''),
            CDS('CDS:ENSP00000504602_4', 41471806, 41471981, 1, ''),
            CDS('CDS:ENSP00000504602_5', 41473325, 41473496, 2, ''),
            CDS('CDS:ENSP00000504602_6', 41476577, 41476620, 1, ''),
            CDS('CDS:ENSP00000504602_7', 41479172, 41479282, 1, ''),
            CDS('CDS:ENSP00000504602_8', 41480476, 41480602, 2, ''),
            CDS('CDS:ENSP00000504602_9', 41488394, 41488513, 2, ''),
            CDS('CDS:ENSP00000504602_10', 41489507, 41489593, 2, ''),
            CDS('CDS:ENSP00000504602_11', 41494356, 41494578, 0, ''),
            CDS('CDS:ENSP00000504602_12', 41498119, 41498189, 2, ''),
            CDS('CDS:ENSP00000504602_13', 41511863, 41511965, 0, '')
        ]
        self.assertEqual(result, cds)


    def test_get_parent_transcript(self):
        """Test case for extracting parent transcript"""
        feature_id = "exon_1012058"
        result = self.db.get_parent_transcript(feature_id)
        transcript = Transcript(
            transcript_id='transcript:ENST00000679263',
            transcript_biotype='protein_coding',
            gene_name='TMPRSS2',
            gene_biotype='protein_coding',
            description='transmembrane serine protease 2 [Source:HGNC Symbol Acc:HGNC:11876]'
        )
        self.assertEqual(result, transcript)
