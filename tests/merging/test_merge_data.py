from bin.merge_data import FusionSummary

class TestFusionSummary:
    input_fusions = "tests/merging/Detected_Fusions.csv"
    input_fusion_context_seqs = "tests/merging/annotated_fusions.csv"
    input_requant_counts = "tests/merging/quantification.tsv"
    input_read_stats = "tests/merging/Log.final.out"
    output_table = "tests/merging/fusions.csv"
    fusion_tools = "fusioncatcher,starfusion,arriba"
    fs = FusionSummary(
        input_fusions,
        input_fusion_context_seqs,
        input_requant_counts,
        input_read_stats,
        output_table,
        fusion_tools
    )


    def test_normalize_counts_cpm(self):
        assert self.fs.normalize_counts_cpm(100) == 0


    def test_load_context_seqs(self):
        assert self.fs.load_context_seqs() == ""


    def test_load_detected_fusions(self):
        assert self.fs.load_detected_fusions() == ""


    def test_load_requant_counts(self):
        assert self.fs.load_requant_counts() == ""
