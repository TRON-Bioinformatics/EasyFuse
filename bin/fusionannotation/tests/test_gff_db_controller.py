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
        # Minimal test DB:
        # GTF file created with `ShrinkGenome` package
        # GFF file generated from GTF file with `AGAT` package
        self.db = GffDbController(
            "tests/fusionannotation/input_data/Homo_sapiens.GRCh38.110_miniannotation.gff3.db"
        )


    def test_get_features_overlapping_position_exon_negative_strand(self):
        """Test case for extracting features overlapping the breakpoint on the negative strand."""
        # 21:41494380:-_7:13935843:- -> TMPRSS2_ETV1
        bp = Breakpoint("21", 30157, "-")
        result = self.db.get_features_overlapping_position(bp, "exon")
        exons = [
            Exon(exon_id='ENSE00003788834', start=30157, stop=30379, transcript_id='ENST00000679263'),
            Exon(exon_id='agat-exon-1358', start=30157, stop=30379, transcript_id='ENST00000679181'),
            Exon(exon_id='agat-exon-1368', start=30157, stop=30379, transcript_id='ENST00000676973'),
            Exon(exon_id='agat-exon-1380', start=30157, stop=30379, transcript_id='ENST00000679054'),
            Exon(exon_id='agat-exon-1393', start=30157, stop=30379, transcript_id='ENST00000678348'),
            Exon(exon_id='agat-exon-1406', start=30157, stop=30379, transcript_id='ENST00000332149'),
            Exon(exon_id='agat-exon-1443', start=30157, stop=30379, transcript_id='ENST00000677680'),
            Exon(exon_id='agat-exon-1450', start=30157, stop=30379, transcript_id='ENST00000678171'),
            Exon(exon_id='agat-exon-1420', start=30157, stop=30379, transcript_id='ENST00000679016'),
            Exon(exon_id='agat-exon-1471', start=30157, stop=30379, transcript_id='ENST00000398585'),
            Exon(exon_id='agat-exon-1483', start=30157, stop=30379, transcript_id='ENST00000458356'),
            Exon(exon_id='ENSE00003909070', start=30157, stop=30379, transcript_id='ENST00000678617'),
            Exon(exon_id='agat-exon-1504', start=30157, stop=30379, transcript_id='ENST00000424093'),
            Exon(exon_id='agat-exon-1515', start=30157, stop=30379, transcript_id='ENST00000454499'),
            Exon(exon_id='agat-exon-1529', start=30157, stop=30379, transcript_id='ENST00000455813')
        ]
        self.assertEqual(result, exons)


    def test_get_features_overlapping_position_exon_positive_strand(self):
        """Test case for extracting features overlapping the breakpoint on the positive strand."""
        bp = Breakpoint("7", 257882, "+")
        result = self.db.get_features_overlapping_position(bp, "exon")
        exons = [
            Exon(exon_id='ENSE00003936914', start=244655, stop=294171, transcript_id='ENST00000691309'),
            Exon(exon_id='ENSE00001886368', start=247896, stop=274015, transcript_id='ENST00000476279'),
            Exon(exon_id='ENSE00001839092', start=247913, stop=273180, transcript_id='ENST00000461457'),
            Exon(exon_id='agat-exon-167', start=257777, stop=257882, transcript_id='ENST00000680072'),
            Exon(exon_id='agat-exon-118', start=257777, stop=257882, transcript_id='ENST00000680513'),
            Exon(exon_id='ENSE00003608380', start=257777, stop=257882, transcript_id='ENST00000680534'),
            Exon(exon_id='agat-exon-217', start=257777, stop=257882, transcript_id='ENST00000679821'),
            Exon(exon_id='agat-exon-264', start=257777, stop=257882, transcript_id='ENST00000681722'),
            Exon(exon_id='agat-exon-311', start=257777, stop=257882, transcript_id='ENST00000680181'),
            Exon(exon_id='agat-exon-361', start=257777, stop=257882, transcript_id='ENST00000356239'),
            Exon(exon_id='agat-exon-410', start=257777, stop=257882, transcript_id='ENST00000679521'),
            Exon(exon_id='agat-exon-513', start=257777, stop=257882, transcript_id='ENST00000679457'),
            Exon(exon_id='ENSE00003626584', start=257777, stop=257882, transcript_id='ENST00000679474'),
            Exon(exon_id='agat-exon-594', start=257777, stop=257882, transcript_id='ENST00000680952'),
            Exon(exon_id='agat-exon-687', start=257777, stop=257882, transcript_id='ENST00000681412'),
            Exon(exon_id='agat-exon-735', start=257777, stop=257882, transcript_id='ENST00000679448'),
            Exon(exon_id='agat-exon-783', start=257777, stop=257882, transcript_id='ENST00000680766'),
            Exon(exon_id='agat-exon-885', start=257777, stop=257882, transcript_id='ENST00000359028'),
            Exon(exon_id='agat-exon-912', start=257777, stop=257882, transcript_id='ENST00000491695'),
            Exon(exon_id='agat-exon-934', start=257777, stop=257882, transcript_id='ENST00000394534'),
            Exon(exon_id='agat-exon-977', start=257777, stop=257882, transcript_id='ENST00000680365'),
            Exon(exon_id='agat-exon-954', start=257777, stop=257882, transcript_id='ENST00000681216'),
            Exon(exon_id='agat-exon-993', start=257777, stop=257882, transcript_id='ENST00000487692'),
            Exon(exon_id='agat-exon-1001', start=257777, stop=257882, transcript_id='ENST00000487258'),
            Exon(exon_id='agat-exon-1008', start=257777, stop=257882, transcript_id='ENST00000486313'),
            Exon(exon_id='ENSE00003912099', start=257777, stop=259537, transcript_id='ENST00000680047'),
            Exon(exon_id='ENSE00001956879', start=257822, stop=257882, transcript_id='ENST00000463118')
        ]
        self.assertEqual(result, exons)


    def test_get_features_overlapping_position_cds(self):
        """Test case for extracting features overlapping the breakpoint."""
        # Original locus: 21:41494380:- -> TMPRSS2
        bp = Breakpoint("21", 30157, "-")
        result = self.db.get_features_overlapping_position(bp, "CDS")
        cds = [
            CDS(cds_id='agat-cds-1367', start=30157, stop=30379, frame=0, transcript_id='ENST00000679263'),
            CDS(cds_id='agat-cds-1380', start=30157, stop=30379, frame=0, transcript_id='ENST00000679181'),
            CDS(cds_id='agat-cds-1391', start=30157, stop=30379, frame=0, transcript_id='ENST00000676973'),
            CDS(cds_id='agat-cds-1404', start=30157, stop=30379, frame=0, transcript_id='ENST00000679054'),
            CDS(cds_id='agat-cds-1417', start=30157, stop=30379, frame=0, transcript_id='ENST00000678348'),
            CDS(cds_id='agat-cds-1430', start=30157, stop=30379, frame=0, transcript_id='ENST00000332149'),
            CDS(cds_id='agat-cds-1465', start=30157, stop=30379, frame=0, transcript_id='ENST00000677680'),
            CDS(cds_id='agat-cds-1472', start=30157, stop=30379, frame=0, transcript_id='ENST00000678171'),
            CDS(cds_id='agat-cds-1443', start=30157, stop=30379, frame=0, transcript_id='ENST00000679016'),
            CDS(cds_id='agat-cds-1489', start=30157, stop=30379, frame=0, transcript_id='ENST00000398585'),
            CDS(cds_id='agat-cds-1502', start=30157, stop=30379, frame=0, transcript_id='ENST00000458356'),
            CDS(cds_id='agat-cds-1515', start=30157, stop=30379, frame=0, transcript_id='ENST00000424093'),
            CDS(cds_id='agat-cds-1527', start=30157, stop=30379, frame=0, transcript_id='ENST00000454499'),
            CDS(cds_id='agat-cds-1540', start=30157, stop=30379, frame=0, transcript_id='ENST00000455813')
        ]
        self.assertEqual(result, cds)


    def test_get_exons_from_transcript(self):
        """Test case for extracting exons from transcript."""
        transcript_id = "ENST00000679263"
        result = self.db.get_features_from_transcript(transcript_id, "exon")
        exons = [
            Exon(exon_id='ENSE00003907476', start=101, stop=1954, transcript_id=''),
            Exon(exon_id='ENSE00003786558', start=3535, stop=3687, transcript_id=''),
            Exon(exon_id='ENSE00001324661', start=4197, stop=4339, transcript_id=''),
            Exon(exon_id='ENSE00001309041', start=6449, stop=6544, transcript_id=''),
            Exon(exon_id='ENSE00001310536', start=7607, stop=7782, transcript_id=''),
            Exon(exon_id='ENSE00001291248', start=9126, stop=9297, transcript_id=''),
            Exon(exon_id='ENSE00001319118', start=12378, stop=12421, transcript_id=''),
            Exon(exon_id='ENSE00001296879', start=14973, stop=15083, transcript_id=''),
            Exon(exon_id='ENSE00001328752', start=16277, stop=16403, transcript_id=''),
            Exon(exon_id='ENSE00001308618', start=24195, stop=24314, transcript_id=''),
            Exon(exon_id='ENSE00003500399', start=25308, stop=25394, transcript_id=''),
            Exon(exon_id='ENSE00003788834', start=30157, stop=30379, transcript_id=''),
            Exon(exon_id='ENSE00003523611', start=33920, stop=33990, transcript_id=''),
            Exon(exon_id='ENSE00003905802', start=47664, stop=47909, transcript_id='')
        ]
        self.assertEqual(result, exons)


    def test_get_cds_from_transcripts(self):
        """Test case for extracting cds from transcript."""
        transcript_id = "ENST00000679263"
        result = self.db.get_features_from_transcript(transcript_id, "CDS")
        cds = [
            CDS(cds_id='agat-cds-1378', start=1943, stop=1954, frame=0, transcript_id=''),
            CDS(cds_id='agat-cds-1377', start=3535, stop=3687, frame=0, transcript_id=''),
            CDS(cds_id='agat-cds-1376', start=4197, stop=4339, frame=2, transcript_id=''),
            CDS(cds_id='agat-cds-1375', start=6449, stop=6544, frame=2, transcript_id=''),
            CDS(cds_id='agat-cds-1374', start=7607, stop=7782, frame=1, transcript_id=''),
            CDS(cds_id='agat-cds-1373', start=9126, stop=9297, frame=2, transcript_id=''),
            CDS(cds_id='agat-cds-1372', start=12378, stop=12421, frame=1, transcript_id=''),
            CDS(cds_id='agat-cds-1371', start=14973, stop=15083, frame=1, transcript_id=''),
            CDS(cds_id='agat-cds-1370', start=16277, stop=16403, frame=2, transcript_id=''),
            CDS(cds_id='agat-cds-1369', start=24195, stop=24314, frame=2, transcript_id=''),
            CDS(cds_id='agat-cds-1368', start=25308, stop=25394, frame=2, transcript_id=''),
            CDS(cds_id='agat-cds-1367', start=30157, stop=30379, frame=0, transcript_id=''),
            CDS(cds_id='agat-cds-1366', start=33920, stop=33990, frame=2, transcript_id=''),
            CDS(cds_id='agat-cds-1365', start=47664, stop=47766, frame=0, transcript_id='')
        ]
        self.assertEqual(result, cds)


    def test_get_parent_transcript(self):
        """Test case for extracting parent transcript."""
        feature_id = "ENSE00003907476"
        result = self.db.get_parent_transcript(feature_id)
        transcript = Transcript(
            transcript_id='ENST00000679263',
            transcript_biotype='protein_coding',
            gene_name='TMPRSS2',
            gene_biotype='protein_coding',
            description='NA'
        )
        self.assertEqual(result, transcript)
