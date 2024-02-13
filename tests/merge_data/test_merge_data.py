from bin.merge_data import FusionSummary

class TestFusionSummary:
    input_fusions = "tests/merge_data/Detected_Fusions.csv"
    input_fusion_context_seqs = "tests/merge_data/annotated_fusions.csv"
    input_requant_counts = "tests/merge_data/quantification.tsv"
    input_read_stats = "tests/merge_data/Log.final.out"
    output_table = "tests/merge_data/fusions.csv"
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
        #print(self.fs.normalize_counts_cpm(100))
        assert round(float(self.fs.normalize_counts_cpm(100))) == 4946


    def test_load_context_seqs(self):
        #print(self.fs.load_context_seqs())
        assert self.fs.load_context_seqs()['ETV6_12:11869969:+_ENST00000396373_NTRK3_15:87940753:-_ENST00000629765'] == {
            'BPID': '12:11869969:+_15:87940753:-', 
            'Fusion_Gene': 'ETV6_NTRK3', 
            'Breakpoint1': '12:11869969:+', 
            'Breakpoint2': '15:87940753:-', 
            'FTID': 'ETV6_12:11869969:+_ENST00000396373_NTRK3_15:87940753:-_ENST00000629765', 
            'context_sequence_id': 'b21d222ed947abeb', 
            'context_sequence_100_id': '158ed13525ea27d6', 
            'type': 'trans_inv', 
            'exon_nr': '11', 
            'exon_starts': '11869424*87940623', 
            'exon_ends': '11869969*87940753', 
            'exon_boundary1': 'right_boundary', 
            'exon_boundary2': 'right_boundary', 
            'exon_boundary': 'both', 
            'bp1_frame': '2', 
            'bp2_frame': '2', 
            'frame': 'in_frame', 
            'context_sequence': 'CCCCTGGACAACATGATCCGCCGCCTCTCCCCGGCTGAGAGAGCTCAGGGACCCAGGCCGCACCAGGAGAACAACCACCAGGAGTCCTACCCTCTGTCAGTGTCTCCCATGGAGAATAATCACTGCCCAGCGTCCTCCGAGTCCCACCCGAAGCCATCCAGCCCCCGGCAGGAGAGCACACGCGTGATCCAGCTGATGCCCAGCCCCATCATGCACCCTCTGATCCTGAACCCCCGGCACTCCGTGGATTTCAAACAGTCCAGGCTCTCCGAGGACGGGCTGCATAGGGAAGGGAAGCCCATCAACCTCTCTCATCGGGAAGACCTGGCTTACATGAACCACATCATGGTCTCTGTCTCCCCGCCTGAAGAGCACGCCATGCCCATTGGGAGAATAGCAGATGTGCAGCACATTAAGAGGAGAGACATCGTGCTGAAGCGAGAACTGGGTGAGGGAGCCTTTGGAAAGGTCTTCCTGGCCGAGTGCTACAACCTCAGCCCGACCAAGGACAAGATGCTTGTGGCTGTGAAGGCCCTGAAGGATCCCACCCTGGCTGCCCGGAAGGATTTCCAGAGGGAGGCCGAGCTGCTCACCAACCTGCAGCATGAGCACATTGTCAAGTTCTATGGAGTGTGCGGCGATGGGGACCCCCTCATCATGGTCTTTGAATACATGAAGCATGGAGACCTGAATAAGTTCCTCAGGGCCCATGGGCCAGATGCAATGATCCTTGTGGATGGACAGCCACGCCAGGCCAAGGGTGAGCTGGGGCTCTCCCAAATGCTCCACATTGCCAGTCA', 
            'context_sequence_bp': '400', 
            'neo_peptide_sequence': 'PPEEHAMPIGRIADVQHIKRRDIVLK', 
            'neo_peptide_sequence_bp': '13.3'
        }


    def test_load_detected_fusions(self):
        assert self.fs.load_detected_fusions() == {
            '21:41494380:-_7:13935843:-': {
                'fusioncatcher': {'junc': '16', 'span': '367'}
            }, 
            '22:29287134:+_12:50814280:+': {
                'fusioncatcher': {'junc': '23', 'span': '363'}, 
                'starfusion': {'junc': '129', 'span': '367'}, 
                'arriba': {'junc': '49', 'span': '299'}
            }, 
            '7:92003235:+_7:140787584:-': {
                'fusioncatcher': {'junc': '21', 'span': '340'}, 
                'starfusion': {'junc': '120', 'span': '348'}, 
                'arriba': {'junc': '113', 'span': '299'}
            }, 
            '12:11869969:+_15:87940753:-': {
                'fusioncatcher': {'junc': '23', 'span': '299'}, 
                'starfusion': {'junc': '160', 'span': '390'}, 
                'arriba': {'junc': '130', 'span': '300'}
            }, 
            '8:42968214:+_10:43116584:+': {
                'fusioncatcher': {'junc': '20', 'span': '298'}, 
                'starfusion': {'junc': '94', 'span': '276'}, 
                'arriba': {'junc': '88', 'span': '275'}
            }, 
            '22:29287134:+_11:128807180:+': {
                'fusioncatcher': {'junc': '22', 'span': '282'}, 
                'starfusion': {'junc': '140', 'span': '333'}, 
                'arriba': {'junc': '70', 'span': '299'}
            }, 
            '19:15254152:-_15:34347969:+': {
                'fusioncatcher': {'junc': '20', 'span': '244'}, 
                'starfusion': {'junc': '132', 'span': '260'}, 
                'arriba': {'junc': '82', 'span': '232'}
            }, 
            '5:150404680:-_6:117324415:-': {
                'fusioncatcher': {'junc': '22', 'span': '152'}, 
                'starfusion': {'junc': '125', 'span': '140'}, 
                'arriba': {'junc': '75', 'span': '139'}
            }, 
            '21:41494375:-_7:13935838:-': {
                'starfusion': {'junc': '76', 'span': '321'}, 
                'arriba': {'junc': '59', 'span': '300'}
            }, 
            '2:42301394:+_2:29223584:-': {
                'starfusion': {'junc': '61', 'span': '256'}, 
                'arriba': {'junc': '52', 'span': '221'}
            }
        }


    def test_load_requant_counts(self):
        #print(self.fs.load_requant_counts()['b21d222ed947abeb'])
        assert self.fs.load_requant_counts()['b21d222ed947abeb'] == {
            'ft': {
                'name': 'b21d222ed947abeb_400_ft', 
                'pos': '400', 
                'junc': '205', 
                'span': '296', 
                'anch': '25', 
                'a': '521', 
                'b': '378'
            }, 
            'wt1': {
                'name': 'b21d222ed947abeb_400_wt1', 
                'pos': '400', 
                'junc': '0', 
                'span': '0', 
                'anch': '0', 
                'a': '2', 
                'b': '0'
            }, 
            'wt2': {
                'name': 'b21d222ed947abeb_400_wt2', 
                'pos': '400', 
                'junc': '0', 
                'span': '0', 
                'anch': '0', 
                'a': '0', 
                'b': '6'
            }
        }
