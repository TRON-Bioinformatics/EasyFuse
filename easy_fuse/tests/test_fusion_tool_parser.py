import unittest

import easy_fuse.tests
import pkg_resources

from easy_fuse.fusiontoolparser_helper import parse_fusioncatcher_results
from easy_fuse.fusiontoolparser_helper import parse_starfusion_results
from easy_fuse.fusiontoolparser_helper import parse_mapsplice_results
from easy_fuse.fusiontoolparser_helper import parse_infusion_results
from easy_fuse.fusiontoolparser_helper import parse_soapfuse_results


class TestFusionToolParser(unittest.TestCase):
    def setUp(self):
        self.fusion_catcher1 = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/fusioncatcher_1.txt"
        )
        self.fusion_catcher2 = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/fusioncatcher_2.txt"
        )
        self.star_fusion = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/starfusion.txt"
        )
        self.mapsplice = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/mapsplice.txt"
        )
        self.infusion = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/infusion.txt"
        )
        self.soapfuse = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/soapfuse.txt"
        )

    def test_parse_fusioncatcher_results(self):
        result_dict = {
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "23",
                "340",
            ],
            "21:41494380:-_7:13935843:-": [
                "TMPRSS2_ETV1",
                "21:41494380:-",
                "7:13935843:-",
                "16",
                "331",
            ],
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "21",
                "308",
            ],
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "23",
                "272",
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "22",
                "268",
            ],
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "20",
                "260",
            ],
            "19:15254152:-_15:34347969:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347969:+",
                "20",
                "213",
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_ROS1",
                "5:150404680:-",
                "6:117324415:-",
                "22",
                "124",
            ],
        }
        self.assertEqual(
            parse_fusioncatcher_results(self.fusion_catcher1, self.fusion_catcher2),
            result_dict,
        )

    def test_parse_starfusion_results(self):
        result_dict = {
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "130",
                "393",
            ],
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "97",
                "288",
            ],
            "21:41494375:-_7:13935838:-": [
                "TMPRSS2_ETV1",
                "21:41494375:-",
                "7:13935838:-",
                "77",
                "364",
            ],
            "19:15254152:-_15:34347969:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347969:+",
                "99",
                "260",
            ],
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "91",
                "277",
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "99",
                "236",
            ],
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "55",
                "197",
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_RP1-179P9.3",
                "5:150404680:-",
                "6:117324415:-",
                "34",
                "67",
            ],
        }
        self.assertEqual(parse_starfusion_results(self.star_fusion), result_dict)

    def test_parse_mapsplice_results(self):
        result_dict = {
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "139",
                "468",
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "153",
                "476",
            ],
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "167",
                "618",
            ],
            "19:15254152:-_15:34347969:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347969:+",
                "182",
                "414",
            ],
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "225",
                "566",
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_ROS1",
                "5:150404680:-",
                "6:117324415:-",
                "179",
                "168",
            ],
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "165",
                "606",
            ],
        }
        self.assertEqual(parse_mapsplice_results(self.mapsplice), result_dict)

    def test_parse_infusion_results(self):
        result_dict = {
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "290",
                "0",
            ],
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "254",
                "0",
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "230",
                "0",
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_ROS1",
                "5:150404680:-",
                "6:117324415:-",
                "201",
                "0",
            ],
            "19:15254152:-_15:34347964:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347964:+",
                "186",
                "0",
            ],
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "190",
                "0",
            ],
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "151",
                "0",
            ],
            "21:41494375:-_7:13935844:-": [
                "TMPRSS2_ETV1",
                "21:41494375:-",
                "7:13935844:-",
                "125",
                "0",
            ],
        }

        self.assertEqual(parse_infusion_results(self.infusion), result_dict)

    def test_parse_soapfuse_results(self):

        result_dict = {
            "7:92003235:+_7:140787584:-": [
                "AKAP9_BRAF",
                "7:92003235:+",
                "7:140787584:-",
                "281",
                "364",
            ],
            "19:15254152:-_15:34347969:+": [
                "BRD4_NUTM1",
                "19:15254152:-",
                "15:34347969:+",
                "153",
                "229",
            ],
            "5:150404680:-_6:117324415:-": [
                "CD74_ROS1",
                "5:150404680:-",
                "6:117324415:-",
                "193",
                "138",
            ],
            "12:11869969:+_15:87940753:-": [
                "ETV6_NTRK3",
                "12:11869969:+",
                "15:87940753:-",
                "245",
                "321",
            ],
            "22:29287134:+_12:50814280:+": [
                "EWSR1_ATF1",
                "22:29287134:+",
                "12:50814280:+",
                "156",
                "329",
            ],
            "22:29287134:+_11:128807180:+": [
                "EWSR1_FLI1",
                "22:29287134:+",
                "11:128807180:+",
                "188",
                "301",
            ],
            "8:42968214:+_10:43116584:+": [
                "HOOK3_RET",
                "8:42968214:+",
                "10:43116584:+",
                "118",
                "254",
            ],
            "21:41494381:-_7:13935844:-": [
                "TMPRSS2_ETV1",
                "21:41494381:-",
                "7:13935844:-",
                "133",
                "323",
            ],
        }

        self.assertEqual(parse_soapfuse_results(self.soapfuse), result_dict)


if __name__ == "__main__":
    unittest.main()
