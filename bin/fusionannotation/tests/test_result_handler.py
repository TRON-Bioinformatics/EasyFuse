"""
Tests for result handler module.
"""

import unittest

# pylint: disable=E0401
from Bio.Seq import Seq # type: ignore

from fusionannotation.src.cds import CDS
from fusionannotation.src.exon import Exon
from fusionannotation.src.transcript import Transcript
from fusionannotation.src.breakpoint import Breakpoint
from fusionannotation.src.fusion_transcript import FusionTranscript
from fusionannotation.src.result_handler import ResultHandler


class TestResultHandler(unittest.TestCase):
    """
    Provides unit tests for the sequence validation module.
    """

    def setUp(self):
        cds_seqs = {
            '21': [
                {
                    Exon(exon_id='exon_1', start=30157, stop=30179, transcript_id='ENST000A'),
                    CDS(cds_id='cds_1', start=30157, stop=30179, frame=0, transcript_id='ENST000A')
                },
                [
                    Seq('TTGAGGTGCACACTGTCCCGGAT'),
                    Seq('TTGAGGTGCACACTGTCCCGGAT')
                ]
            ],
            '7': [
                {
                    Exon(exon_id='exon_2', start=248000, stop=248030, transcript_id='ENST000B'),
                    CDS(cds_id='cds_2', start=248000, stop=248030, frame=0, transcript_id='ENST000B')
                },
                [
                    Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT'),
                    Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT')
                ]
            ]
        }
        self.result_handler = ResultHandler(cds_seqs, 100000)


    def test_generate_temp_dict(self):
        """Test case for generation of temporary dictionary"""
        exon = Exon(exon_id='exon_1', start=30157, stop=30179, transcript_id='ENST000A')
        cds = CDS(cds_id='cds_1', start=30157, stop=30179, frame=0, transcript_id='ENST000A')
        tmp_dict = {
            exon: Seq('TTGAGGTGCACACTGTCCCGGAT'),
            cds: Seq('TTGAGGTGCACACTGTCCCGGAT')
        }
        result = self.result_handler.generate_temp_dict('21')
        self.assertEqual(result, tmp_dict)


    def test_generate_fusion_transcript_values(self):
        """Test case to generate result row"""
        wt1 = Transcript("ENST000A", "", "GeneA", "", "")
        wt2 = Transcript("ENST000B", "", "GeneB", "", "")
        bp1 = Breakpoint("21", 30157, "-")
        bp2 = Breakpoint("7", 248000, "+")
        bp1.exon_boundary = "left_boundary"
        bp2.exon_boundary = "right_boundary"
        bp1.cds_boundary = "left_boundary"
        bp2.cds_boundary = "right_boundary"
        wt1.frame_at_start = 0
        wt2.frame_at_start = 0
        wt1.exons = [
            Exon(exon_id='exon_1', start=30157, stop=30179, transcript_id='ENST000A'),
        ]
        wt1.cds = [
            CDS(cds_id='cds_1', start=30157, stop=30179, frame=0, transcript_id='ENST000A')
        ]
        wt2.exons = [
            Exon(exon_id='exon_2', start=248000, stop=248030, transcript_id='ENST000B')
        ]
        wt2.cds = [
            CDS(cds_id='cds_2', start=248000, stop=248030, frame=0, transcript_id='ENST000B')
        ]
        wt1.is_good_transcript = set()
        wt2.is_good_transcript = set()
        fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
        fusion_transcript.frame = "in_frame"
        fusion_transcript.fusion_type = "trans"
        actual_result = self.result_handler.generate_fusion_transcript_values(fusion_transcript)
        expected_result = {
            'BPID': '21:30157:-_7:248000:+',
            'Fusion_Gene': 'GeneA_GeneB',
            'Breakpoint1': '21:30157:-',
            'Breakpoint2': '7:248000:+',
            'FTID': 'GeneA_21:30157:-_ENST000A_GeneB_7:248000:+_ENST000B',
            'context_sequence_id': '6fc9ef059658b069',
            'context_sequence_100_id': '6fc9ef059658b069',
            'type': 'trans',
            'exon_nr': '2',
            'ft1_exon_nr': '1',
            'ft2_exon_nr': '0',
            'exon_starts': '30157*248000',
            'exon_ends': '30179*248030',
            'exon_boundary1': 'left_boundary',
            'exon_boundary2': 'right_boundary',
            'exon_boundary': 'no_match',
            'bp1_frame': 'None',
            'bp2_frame': 'None',
            'frame': 'in_frame',
            'context_sequence': Seq('ATCCGGGACAGTGTGCACCTCAATCTATCAAAATAAAAAATAGTGACAGCAAGT'),
            'context_sequence_bp': 1,
            'neo_peptide_sequence': Seq('IRDSVHLNLSK'),
            'neo_peptide_sequence_bp': 0.3,
            'fusion_protein_sequence': Seq('IRDSVHLNLSK'),
            'fusion_protein_sequence_bp': 0.3,
            'context_sequence_wt1': Seq('ATCCGGGACAGTGTGCACCTCAA'),
            'context_sequence_wt2': Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT'),
            'context_sequence_wt1_bp': 1,
            'context_sequence_wt2_bp': 0,
            'context_sequence_100': Seq('ATCCGGGACAGTGTGCACCTCAATCTATCAAAATAAAAAATAGTGACAGCAAGT'),
            'bp1_chr': '21',
            'bp1_pos': 30157,
            'bp1_strand': '-',
            'bp2_chr': '7',
            'bp2_pos': 248000,
            'bp2_strand': '+',
            'cds_boundary1': 'left_boundary',
            'cds_boundary2': 'right_boundary',
            'wt1_exon_pos': [
                Exon(exon_id='exon_1', start=30157, stop=30179, transcript_id='ENST000A')
            ],
            'wt2_exon_pos': [
                Exon(exon_id='exon_2', start=248000, stop=248030, transcript_id='ENST000B')
            ],
            'ft1_exon_pos': [
                Exon(exon_id='exon_1', start=30157, stop=30179, transcript_id='ENST000A')
            ],
            'ft2_exon_pos': [
                Exon(exon_id='exon_2', start=248000, stop=248030, transcript_id='ENST000B')
            ],
            'wt1_cds_pos': [
                CDS(cds_id='cds_1', start=30157, stop=30179, frame=0, transcript_id='ENST000A')
            ],
            'wt2_cds_pos': [
                CDS(cds_id='cds_2', start=248000, stop=248030, frame=0, transcript_id='ENST000B')
            ],
            'ft1_cds_pos': [
                CDS(cds_id='cds_1', start=30157, stop=30179, frame=0, transcript_id='ENST000A')
            ],
            'ft2_cds_pos': [
                CDS(cds_id='cds_2', start=248000, stop=248030, frame=0, transcript_id='ENST000B')
            ],
            'wt1_exon_seqs': [
                Seq('TTGAGGTGCACACTGTCCCGGAT')
            ],
            'wt2_exon_seqs': [
                Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT')
            ],
            'ft1_exon_seqs': [
                Seq('TTGAGGTGCACACTGTCCCGGAT')
            ],
            'ft2_exon_seqs': [
                Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT')
            ],
            'wt1_cds_seqs': [
                Seq('TTGAGGTGCACACTGTCCCGGAT')
            ],
            'wt2_cds_seqs': [
                Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT')
            ],
            'ft1_cds_seqs': [
                Seq('TTGAGGTGCACACTGTCCCGGAT')
            ],
            'ft2_cds_seqs': [
                Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT')
            ],
            'wt1_exon_transcripts': Seq('TTGAGGTGCACACTGTCCCGGAT'),
            'wt2_exon_transcripts': Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT'),
            'ft1_exon_transcripts': Seq('TTGAGGTGCACACTGTCCCGGAT'),
            'ft2_exon_transcripts': Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT'),
            'wt1_cds_transcripts': 'TTGAGGTGCACACTGTCCCGGAT',
            'wt2_cds_transcripts': 'TCTATCAAAATAAAAAATAGTGACAGCAAGT',
            'ft1_cds_transcripts': Seq('TTGAGGTGCACACTGTCCCGGAT'),
            'ft2_cds_transcripts': Seq('TCTATCAAAATAAAAAATAGTGACAGCAAGT'),
            'wt1_peptide': Seq('IRDSVHL'), 'wt2_peptide': Seq('SIKIKNSDSK'),
            'fusion_transcript': Seq('ATCCGGGACAGTGTGCACCTCAATCTATCAAAATAAAAAATAGTGACAGCAAGT'),
            'fusion_peptide': Seq('IRDSVHLNLSK'),
            'wt1_is_good_transcript': {'wt1 seq % 3 != 0'},
            'wt2_is_good_transcript': {'wt2 seq % 3 != 0'},
            'wt1_trans_biotype': '',
            'wt2_trans_biotype': '',
            'wt1_gene_biotype': '',
            'wt2_gene_biotype': '',
            'wt1_description': '',
            'wt2_description': '',
            'wt1_frame_at_start': 0,
            'wt2_frame_at_start': 0,
            'wt1_TSL': None,
            'wt2_TSL': None,
            'wt1_exon_no': 1,
            'wt2_exon_no': 1,
            'ft1_exon_no': 1,
            'ft2_exon_no': 1,
            'wt1_cds_no': 1,
            'wt2_cds_no': 1,
            'ft1_cds_no': 1,
            'ft2_cds_no': 1,
            'wt1_start_stop': '21:30157:30179',
            'wt2_start_stop': '7:248000:248030',
            'annotation_bias': False
        }
        self.assertEqual(actual_result, expected_result)


if __name__ == '__main__':
    unittest.main()
