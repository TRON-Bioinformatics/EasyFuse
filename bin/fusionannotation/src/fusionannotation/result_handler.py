"""
This module contains functions to handle the results of the fusion annotation.
"""

# pylint: disable=E0401
from .exon import get_exon_ranges
from .exon import get_exon_ranges_reverse
from .fusion_transcript import FusionTranscript
from .sequence_handler import get_stranded_seq
from .sequence_handler import get_feature_seqs
from .sequence_handler import concatenate_seqs
from .sequence_handler import get_context_sequence
from .sequence_handler import calc_hash
from .sequence_handler import get_peptide_sequence
from .sequence_handler import get_fusion_transcript_sequence
from .sequence_handler import get_trimmed_seq


class ResultHandler:
    """
    Class to handle the results of the fusion annotation.
    """

    def __init__(self, cds_seqs: dict, context_seq_len: int):
        self.ftid_plus_set = set()
        self.cds_seqs = cds_seqs
        self.context_seq_len = context_seq_len


    def generate_temp_dict(self, chrom: str) -> dict:
        """Generates a temporary dictionary for a chromosome.

        Args:
            cds_seq_dict (dict): Feature dictionary with according sequences
            chrom (str): Chromosome to generate the dictionary for

        Returns:
            dict: Temporary dictionary with cds positions as keys and sequences as values
        """
        return dict(
            zip(
                self.cds_seqs[chrom][0],
                self.cds_seqs[chrom][1],
            )
        )


    def generate_fusion_transcript_values(self, fusion_transcript: FusionTranscript) -> dict:
        """Generates a result row for the fusion transcript.

        Args:
            fusion_transcript (FusionTranscript): Fusion transcript object
            context_seq_len (int): Length of the context sequence

        Returns:
            dict: Result row for the fusion transcript
        """
        wt1 = fusion_transcript.transcript_1
        wt2 = fusion_transcript.transcript_2
        bp1 = fusion_transcript.bp1
        bp2 = fusion_transcript.bp2
        tmp_dict_1 = self.generate_temp_dict(bp1.chrom)
        tmp_dict_2 = self.generate_temp_dict(bp2.chrom)

        wt1_cds_seqs = get_feature_seqs(tmp_dict_1, wt1.cds)
        wt2_cds_seqs = get_feature_seqs(tmp_dict_2, wt2.cds)
        wt1_cds_transcripts = concatenate_seqs(wt1_cds_seqs)
        wt2_cds_transcripts = concatenate_seqs(wt2_cds_seqs)

        ft1_cds_seqs = get_feature_seqs(tmp_dict_1, fusion_transcript.cds_transcript_1)
        ft2_cds_seqs = get_feature_seqs(tmp_dict_2, fusion_transcript.cds_transcript_2)
        ft1_cds_transcripts = concatenate_seqs(ft1_cds_seqs)
        ft2_cds_transcripts = concatenate_seqs(ft2_cds_seqs)

        ft1_exon_seqs = get_feature_seqs(tmp_dict_1, fusion_transcript.exons_transcript_1)
        ft2_exon_seqs = get_feature_seqs(tmp_dict_2, fusion_transcript.exons_transcript_2)
        ft1_exon_transcripts = concatenate_seqs(ft1_exon_seqs)
        ft2_exon_transcripts = concatenate_seqs(ft2_exon_seqs)

        wt1_exon_seqs = get_feature_seqs(tmp_dict_1, wt1.exons)
        wt2_exon_seqs = get_feature_seqs(tmp_dict_2, wt2.exons)
        wt1_exon_transcripts = concatenate_seqs(wt1_exon_seqs)
        wt2_exon_transcripts = concatenate_seqs(wt2_exon_seqs)

        wt1_peptide = get_peptide_sequence(wt1_cds_transcripts, bp1, wt1.frame_at_start)
        wt2_peptide = get_peptide_sequence(wt2_cds_transcripts, bp2, wt2.frame_at_start)

        fusion_transcript_seq = get_fusion_transcript_sequence(
            ft1_cds_transcripts,
            ft2_cds_transcripts,
            ft2_exon_transcripts,
            bp1.strand,
            bp2.strand,
            fusion_transcript.frame
        )

        # translate and trim the fusion sequence around the breakpoint
        fusion_peptide = fusion_transcript_seq[wt1.frame_at_start:].translate(
            table=1, to_stop=True
        )

        # the fusion starts in the fusion transcript with the end of the first part
        bp_in_fusion_nt_ex = len(fusion_transcript.exons_transcript_2)

        context_sequence_stranded = get_context_sequence(
            ft1_exon_transcripts,
            ft2_exon_transcripts,
            bp1.strand,
            bp2.strand
        )
        context_sequence = get_trimmed_seq(
            context_sequence_stranded,
            max(0, bp_in_fusion_nt_ex - self.context_seq_len),
            (bp_in_fusion_nt_ex + self.context_seq_len)
        )

        context_sequence_100 = get_trimmed_seq(
            context_sequence,
            max(0, bp_in_fusion_nt_ex - 100),
            (bp_in_fusion_nt_ex + 100)
        )

        ctx_id = calc_hash(context_sequence)
        context_sequence_100_id = calc_hash(context_sequence_100)

        # the breakpoint in wt1 must be identical to the breakpoint in the fusion
        bp_in_wt1_nt_ex = bp_in_fusion_nt_ex
        context_sequence_wt1_stranded = get_stranded_seq(
            wt1_exon_transcripts,
            bp1.strand
        )
        context_sequence_wt1 = get_trimmed_seq(
            context_sequence_wt1_stranded,
            max(0, bp_in_wt1_nt_ex - self.context_seq_len),
            (bp_in_wt1_nt_ex + self.context_seq_len)
        )

        # the breakpoint in wt2 is at the position where the ft2 transcripts start
        bp_in_wt2_nt_ex = len(wt2_exon_transcripts) - len(ft2_exon_transcripts)
        context_sequence_wt2_stranded = get_stranded_seq(
            wt2_exon_transcripts,
            bp2.strand
        )
        context_sequence_wt2 = get_trimmed_seq(
            context_sequence_wt2_stranded,
            max(0, bp_in_wt2_nt_ex - self.context_seq_len),
            bp_in_wt2_nt_ex + self.context_seq_len
        )
        context_sequence_bp = min(self.context_seq_len, bp_in_fusion_nt_ex)
        context_sequence_wt1_bp = min(self.context_seq_len, bp_in_wt1_nt_ex)
        context_sequence_wt2_bp = min(self.context_seq_len, bp_in_wt2_nt_ex)

        # bp pos in the fusion
        # the fusion starts in the fusion transcript with the end of the first part
        bp_in_fusion_nt = len(fusion_transcript.cds_transcript_2)
        # the respective bp pos in the aa seq, whereas x.0 / x.3 / x.6
        # indicate that the breakpoint is at the beginning / after first base /
        # after second base of the underlying codon
        bp_in_fusion_aa = (
            (bp_in_fusion_nt - wt1.frame_at_start) * 10.0 // 3
        ) / 10.0
        neo_peptide_sequence = fusion_peptide[
            max(0, int(bp_in_fusion_aa) - 13) :
        ]

        fusion_protein_sequence = fusion_peptide

        fusion_protein_sequence_bp = round(bp_in_fusion_aa, 1)

        if (int(bp_in_fusion_aa) - 13) < 0:
            neo_peptide_sequence_bp = round(bp_in_fusion_aa, 1)
        else:
            neo_peptide_sequence_bp = round(
                bp_in_fusion_aa - (int(bp_in_fusion_aa) - 13), 1
            )
        if fusion_transcript.frame == "in_frame":
            neo_peptide_sequence = neo_peptide_sequence[
                : (int(neo_peptide_sequence_bp) + 13)
            ]


        ftid = fusion_transcript.get_ftid()
        ftid_plus = f"{ftid}_{ctx_id}"
        if ftid_plus in self.ftid_plus_set:
            annotation_bias = True
        else:
            annotation_bias = False
            self.ftid_plus_set.add(ftid_plus)

        # set exon starts/ends
        # urla: needs revision and simplification, but is working...
        neo_pep_nuc_length = len(neo_peptide_sequence) * 3
        neo_pep_until_bp_nuc = int(round(neo_peptide_sequence_bp * 3, 0))

        if bp1.strand == "+":
            (exon_starts_1, exon_ends_1) = get_exon_ranges_reverse(
                fusion_transcript.exons_transcript_1,
                neo_pep_until_bp_nuc
            )
        else:
            (exon_starts_1, exon_ends_1) = get_exon_ranges(
                fusion_transcript.exons_transcript_1,
                neo_pep_until_bp_nuc
            )

        if bp2.strand == "+":
            (exon_starts_2, exon_ends_2) = get_exon_ranges(
                fusion_transcript.exons_transcript_2,
                neo_pep_nuc_length - neo_pep_until_bp_nuc
            )
        else:
            (exon_starts_2, exon_ends_2) = get_exon_ranges_reverse(
                fusion_transcript.exons_transcript_2,
                neo_pep_nuc_length - neo_pep_until_bp_nuc
            )

        exon_starts = f"{exon_starts_1}*{exon_starts_2}"
        exon_ends = f"{exon_ends_1}*{exon_ends_2}"

        # experimental
        wt1_is_good_transcript = wt1.flags
        if len(wt1_cds_transcripts) % 3 != 0:
            wt1_is_good_transcript.add(
                "wt1 seq % 3 != 0"
            )
        wt2_is_good_transcript = wt2.flags
        if len(wt2_cds_transcripts) % 3 != 0:
            wt2_is_good_transcript.add(
                "wt2 seq % 3 != 0"
            )

        wt1_cds_transcripts = str(
            wt1_cds_transcripts
        )
        wt2_cds_transcripts = str(
            wt2_cds_transcripts
        )

        table_row = {
            "BPID": fusion_transcript.get_bpid(),
            "Fusion_Gene": fusion_transcript.get_fusion_gene_name(),
            "Breakpoint1": str(bp1),
            "Breakpoint2": str(bp2),
            "FTID": ftid,
            "context_sequence_id": ctx_id,
            "context_sequence_100_id": context_sequence_100_id,
            "type": fusion_transcript.fusion_type,
            "exon_nr": str(fusion_transcript.get_exon_nr()),
            "ft1_exon_nr": str(len(fusion_transcript.exons_transcript_1)),
            "ft2_exon_nr": str(len(wt1.exons) - len(fusion_transcript.exons_transcript_2)),
            "exon_starts": exon_starts,
            "exon_ends": exon_ends,
            "exon_boundary1": bp1.exon_boundary,
            "exon_boundary2": bp2.exon_boundary,
            "exon_boundary": fusion_transcript.get_combined_boundary(),
            "bp1_frame": str(wt1.frame_at_bp),
            "bp2_frame": str(wt2.frame_at_bp),
            "frame": fusion_transcript.frame,
            "context_sequence": context_sequence,
            "context_sequence_bp": context_sequence_bp,
            "neo_peptide_sequence": neo_peptide_sequence,
            "neo_peptide_sequence_bp": neo_peptide_sequence_bp,
            "fusion_protein_sequence": fusion_protein_sequence,
            "fusion_protein_sequence_bp": fusion_protein_sequence_bp,
            "context_sequence_wt1": context_sequence_wt1,
            "context_sequence_wt2": context_sequence_wt2,
            "context_sequence_wt1_bp": context_sequence_wt1_bp,
            "context_sequence_wt2_bp": context_sequence_wt2_bp,
            "context_sequence_100": context_sequence_100,
            "bp1_chr": bp1.chrom,
            "bp1_pos": bp1.pos,
            "bp1_strand": bp1.strand,
            "bp2_chr": bp2.chrom,
            "bp2_pos": bp2.pos,
            "bp2_strand": bp2.strand,
            "cds_boundary1": bp1.cds_boundary,
            "cds_boundary2": bp2.cds_boundary,
            "wt1_exon_pos": wt1.exons,
            "wt2_exon_pos": wt2.exons,
            "ft1_exon_pos": fusion_transcript.exons_transcript_1,
            "ft2_exon_pos": fusion_transcript.exons_transcript_2,
            "wt1_cds_pos": wt1.cds,
            "wt2_cds_pos": wt2.cds,
            "ft1_cds_pos": fusion_transcript.cds_transcript_1,
            "ft2_cds_pos": fusion_transcript.cds_transcript_2,
            "wt1_exon_seqs": wt1_exon_seqs,
            "wt2_exon_seqs": wt2_exon_seqs,
            "ft1_exon_seqs": ft1_exon_seqs,
            "ft2_exon_seqs": ft2_exon_seqs,
            "wt1_cds_seqs": wt1_cds_seqs,
            "wt2_cds_seqs": wt2_cds_seqs,
            "ft1_cds_seqs": ft1_cds_seqs,
            "ft2_cds_seqs": ft2_cds_seqs,
            "wt1_exon_transcripts": wt1_exon_transcripts,
            "wt2_exon_transcripts": wt2_exon_transcripts,
            "ft1_exon_transcripts": ft1_exon_transcripts,
            "ft2_exon_transcripts": ft2_exon_transcripts,
            "wt1_cds_transcripts": wt1_cds_transcripts,
            "wt2_cds_transcripts": wt2_cds_transcripts,
            "ft1_cds_transcripts": ft1_cds_transcripts,
            "ft2_cds_transcripts": ft2_cds_transcripts,
            "wt1_peptide": wt1_peptide,
            "wt2_peptide": wt2_peptide,
            "fusion_transcript": fusion_transcript_seq,
            "fusion_peptide": fusion_peptide,
            "wt1_is_good_transcript": wt1_is_good_transcript,
            "wt2_is_good_transcript": wt2_is_good_transcript,
            "wt1_trans_biotype": wt1.transcript_biotype,
            "wt2_trans_biotype": wt2.transcript_biotype,
            "wt1_gene_biotype": wt1.gene_biotype,
            "wt2_gene_biotype": wt2.gene_biotype,
            "wt1_description": wt1.description,
            "wt2_description": wt2.description,
            "wt1_frame_at_start": wt1.frame_at_start,
            "wt2_frame_at_start": wt2.frame_at_start,
            "wt1_TSL": wt1.tsl,
            "wt2_TSL": wt2.tsl,
            "wt1_exon_no": len(wt1.exons),
            "wt2_exon_no": len(wt2.exons),
            "ft1_exon_no": len(fusion_transcript.exons_transcript_1),
            "ft2_exon_no": len(fusion_transcript.exons_transcript_2),
            "wt1_cds_no": len(wt1.cds),
            "wt2_cds_no": len(wt2.cds),
            "ft1_cds_no": len(fusion_transcript.cds_transcript_1),
            "ft2_cds_no": len(fusion_transcript.cds_transcript_2),
            "wt1_start_stop": f"{bp1.chrom}:{wt1.exons[0].start}:{wt1.exons[-1].stop}",
            "wt2_start_stop": f"{bp2.chrom}:{wt2.exons[0].start}:{wt2.exons[-1].stop}",
            "annotation_bias": annotation_bias
        }
        return table_row



    def generate_result_list(self, fusion_transcripts: list) -> list:
        """
        Based on the gathered information, create and add sequences to the table
        urla note: we could have used a pandas df at create functions to
        apply on the columns to do the following,
        but as the following is mainly string manipulation and simple lookups
        should not make a big difference in terms of speed and readability
        we save the last chromosome string to avoid some superfluous re-calculations
        due to annotation biases between different prediction tools,
        it is possible that the same ftid is created twice
        as this is 100% redundant information, we will skip such records.
        In order to guarantee that this is working as intended,
        we introduce the ftid plus which is the ftid plus appended seq hash
        """

        result_list = []
        for fusion_transcript in fusion_transcripts:
            ft_values = self.generate_fusion_transcript_values(fusion_transcript)
            result_list.append(ft_values)
        return result_list
