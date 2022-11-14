import os
import unittest

from fusiontoolparser_helper import parse_fusioncatcher_results
from fusiontoolparser_helper import parse_starfusion_results
from fusiontoolparser_helper import parse_mapsplice_results
from fusiontoolparser_helper import parse_infusion_results
from fusiontoolparser_helper import parse_soapfuse_results

from fusionannotation import FusionAnnotation


TESTDATA_FUSIONCATCHER_1 = os.path.join(os.path.dirname(__file__), "test_case", "results", "fusioncatcher_1.txt")
TESTDATA_FUSIONCATCHER_2 = os.path.join(os.path.dirname(__file__), "test_case", "results", "fusioncatcher_2.txt")
TESTDATA_STARFUSION = os.path.join(os.path.dirname(__file__), "test_case", "results", "starfusion.txt")
TESTDATA_MAPSPLICE = os.path.join(os.path.dirname(__file__), "test_case", "results", "mapsplice.txt")
TESTDATA_INFUSION = os.path.join(os.path.dirname(__file__), "test_case", "results", "infusion.txt")
TESTDATA_SOAPFUSE = os.path.join(os.path.dirname(__file__), "test_case", "results", "soapfuse.txt")

class TestFusionToolParser(unittest.TestCase):

    def test_parse_fusioncatcher_results(self):
        result_dict = {
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "23",
                "340"
            ],
            "21:41494380:-_7:13935843:-": [
                "TMPRSS2_ETV1",
                "21:41494380:-",
                "7:13935843:-",
                "16",
                "331"
            ],
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "21",
                "308"
            ],
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "23",
                "272"
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "22",
                "268"
            ],
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "20",
                "260"
            ],
            "19:15254152:-_15:34347969:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347969:+",
                "20",
                "213"
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_ROS1",
                "5:150404680:-",
                "6:117324415:-",
                "22",
                "124"
            ]
        }
        self.assertEqual(parse_fusioncatcher_results(TESTDATA_FUSIONCATCHER_1, TESTDATA_FUSIONCATCHER_2), result_dict)

    def test_parse_starfusion_results(self):
        result_dict = {
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "130",
                "393"
            ],
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "97",
                "288"
            ],
            "21:41494375:-_7:13935838:-": [
                "TMPRSS2_ETV1",
                "21:41494375:-",
                "7:13935838:-",
                "77",
                "364"
            ],
            "19:15254152:-_15:34347969:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347969:+",
                "99",
                "260"
            ],
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "91",
                "277"
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "99",
                "236"
            ],
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "55",
                "197"
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_RP1-179P9.3",
                "5:150404680:-",
                "6:117324415:-",
                "34",
                "67"
            ]
        }
        self.assertEqual(parse_starfusion_results(TESTDATA_STARFUSION), result_dict)


    def test_parse_mapsplice_results(self):
        result_dict = {
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "139",
                "468"
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "153",
                "476"
            ],
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "167",
                "618"
            ],
            "19:15254152:-_15:34347969:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347969:+",
                "182",
                "414"
            ],
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "225",
                "566"
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_ROS1",
                "5:150404680:-",
                "6:117324415:-",
                "179",
                "168"
            ],
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "165",
                "606"
            ]
        }
        self.assertEqual(parse_mapsplice_results(TESTDATA_MAPSPLICE), result_dict)


    def test_parse_infusion_results(self):
        result_dict = {
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "290",
                "0"
            ],
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "254",
                "0"
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "230",
                "0"
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_ROS1",
                "5:150404680:-",
                "6:117324415:-",
                "201",
                "0"
            ],
            "19:15254152:-_15:34347964:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347964:+",
                "186",
                "0"
            ],
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "190",
                "0"
            ],
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "151",
                "0"
            ],
            "21:41494375:-_7:13935844:-": [
                "TMPRSS2_ETV1",
                "21:41494375:-",
                "7:13935844:-",
                "125",
                "0"
            ]
        }
        
        self.assertEqual(parse_infusion_results(TESTDATA_INFUSION), result_dict)

    def test_parse_soapfuse_results(self):

        result_dict = {
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "281",
                "364"
            ],
            "19:15254152:-_15:34347969:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347969:+",
                "153",
                "229"
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_ROS1",
                "5:150404680:-",
                "6:117324415:-",
                "193",
                "138"
            ],
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "245",
                "321"
            ],
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "156",
                "329"
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "188",
                "301"
            ],
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "118",
                "254"
            ],
            "21:41494381:-_7:13935844:-": [
                "TMPRSS2_ETV1",
                "21:41494381:-",
                "7:13935844:-",
                "133",
                "323"
            ]
        }

        self.assertEqual(parse_soapfuse_results(TESTDATA_SOAPFUSE), result_dict)

        
class TestFusionAnnotation(unittest.TestCase):

    def test_define_type(self):

        self.assertEqual(FusionAnnotation.define_type(500, "1", 200, "+", "2", 600, "-"), "trans_inv")
        self.assertEqual(FusionAnnotation.define_type(500, "1", 200, "+", "1", 800, "+"), "cis_far")
        self.assertEqual(FusionAnnotation.define_type(1000, "1", 200, "+", "1", 800, "+"), "cis_near")
        self.assertEqual(FusionAnnotation.define_type(500, "4", 800, "+", "4", 200, "+"), "cis_trans")
        self.assertEqual(FusionAnnotation.define_type(500, "7", 800, "-", "7", 200, "-"), "cis_far")
        self.assertEqual(FusionAnnotation.define_type(500, "2", 800, "-", "3", 200, "-"), "trans")


    def test_check_exon_overlap(self):
        
        wt1_exon_pos_list = [(1, 100), (200, 300)]
        wt2_exon_pos_list = [(500, 600), (700, 800)]
        self.assertEqual(FusionAnnotation.check_exon_overlap(wt1_exon_pos_list, wt2_exon_pos_list, "cis_near"), False)

        wt1_exon_pos_list = [(1, 100), (200, 300)]
        wt2_exon_pos_list = [(150, 600), (700, 800)]
        self.assertEqual(FusionAnnotation.check_exon_overlap(wt1_exon_pos_list, wt2_exon_pos_list, "cis_near"), True)

        wt1_exon_pos_list = [(50, 100), (200, 300)]
        wt2_exon_pos_list = [(301, 600), (700, 800)]
        self.assertEqual(FusionAnnotation.check_exon_overlap(wt1_exon_pos_list, wt2_exon_pos_list, "cis_far"), False)

        wt1_exon_pos_list = [(400, 600), (800, 1000)]
        wt2_exon_pos_list = [(200, 300), (400, 500)]
        self.assertEqual(FusionAnnotation.check_exon_overlap(wt1_exon_pos_list, wt2_exon_pos_list, "cis_trans"), True)

        wt1_exon_pos_list = [(1, 100), (200, 300)]
        wt2_exon_pos_list = [(500, 600), (700, 800)]
        self.assertEqual(FusionAnnotation.check_exon_overlap(wt1_exon_pos_list, wt2_exon_pos_list, "trans"), False)


if __name__ == "__main__":
    unittest.main()
