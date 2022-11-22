import unittest

from easy_fuse.fusionannotation import FusionAnnotation


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
