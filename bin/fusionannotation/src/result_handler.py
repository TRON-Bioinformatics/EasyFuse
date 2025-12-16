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
# from .sequence_handler import get_fusion_transcript_sequence
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
        """Generates a result row for the fusion transcript."""
        # Full (original) transcripts
        wt1_full = fusion_transcript.transcript_1
        wt2_full = fusion_transcript.transcript_2
        # Truncated (fusion) feature subsets
        bp1 = fusion_transcript.bp1
        bp2 = fusion_transcript.bp2
        ft_exons_1 = fusion_transcript.exons_transcript_1
        ft_exons_2 = fusion_transcript.exons_transcript_2
        ft_cds_1 = fusion_transcript.cds_transcript_1
        ft_cds_2 = fusion_transcript.cds_transcript_2

        # Generate temporary feature dicts for both chromosomes
        tmp_dict_1 = self.generate_temp_dict(bp1.chrom)
        tmp_dict_2 = self.generate_temp_dict(bp2.chrom)

        # Wildtype (full) sequences
        wt1_cds_seqs = get_feature_seqs(tmp_dict_1, wt1_full.cds)
        wt2_cds_seqs = get_feature_seqs(tmp_dict_2, wt2_full.cds)
        wt1_cds_transcripts = concatenate_seqs(wt1_cds_seqs)
        wt2_cds_transcripts = concatenate_seqs(wt2_cds_seqs)

        # Fusion (truncated) CDS sequences
        ft1_cds_seqs = get_feature_seqs(tmp_dict_1, ft_cds_1)
        ft2_cds_seqs = get_feature_seqs(tmp_dict_2, ft_cds_2)
        ft1_cds_transcripts = concatenate_seqs(ft1_cds_seqs)
        # ft2_cds_transcripts does only contain CDS features of the second transcript after the breakpoint
        # In case of out-of-frame fusions, this may not reflect the actual CDS of the fusion transcript
        # as the first CDS may start later in the exon
        ft2_cds_transcripts = concatenate_seqs(ft2_cds_seqs)

        # Fusion (truncated) exon sequences
        ft1_exon_seqs = get_feature_seqs(tmp_dict_1, ft_exons_1)
        ft2_exon_seqs = get_feature_seqs(tmp_dict_2, ft_exons_2)
        ft1_exon_transcripts = concatenate_seqs(ft1_exon_seqs)
        ft2_exon_transcripts = concatenate_seqs(ft2_exon_seqs)

        # Wildtype exon sequences (full)
        wt1_exon_seqs = get_feature_seqs(tmp_dict_1, wt1_full.exons)
        wt2_exon_seqs = get_feature_seqs(tmp_dict_2, wt2_full.exons)
        wt1_exon_transcripts = concatenate_seqs(wt1_exon_seqs)
        wt2_exon_transcripts = concatenate_seqs(wt2_exon_seqs)

        wt1_peptide = get_peptide_sequence(wt1_cds_transcripts, bp1, wt1_full.frame_at_start)
        wt2_peptide = get_peptide_sequence(wt2_cds_transcripts, bp2, wt2_full.frame_at_start)

        # Breakpoint in fusion transcript (exonic)
        bp_in_fusion_nt_ex = len(ft1_exon_transcripts)

        context_sequence_stranded = get_context_sequence(
            ft1_seq = ft1_exon_transcripts,
            ft2_seq = ft2_exon_transcripts,
            bp1_strand = bp1.strand,
            bp2_strand = bp2.strand
        )

        context_sequence = get_trimmed_seq(
            context_sequence_stranded,
            max(0, bp_in_fusion_nt_ex - self.context_seq_len),
            min(len(context_sequence_stranded), bp_in_fusion_nt_ex + self.context_seq_len)
        )

        context_sequence_100 = get_trimmed_seq(
            context_sequence_stranded,
            max(0, bp_in_fusion_nt_ex - 100),
            min(len(context_sequence_stranded), bp_in_fusion_nt_ex + 100)
        )

        ctx_id = calc_hash(context_sequence)
        context_sequence_100_id = calc_hash(context_sequence_100)

        # the breakpoint in wt1 must be identical to the breakpoint in the fusion
        bp_in_wt1_nt_ex = bp_in_fusion_nt_ex
        context_sequence_wt1_stranded = get_stranded_seq(wt1_exon_transcripts, bp1.strand)
        context_sequence_wt1 = get_trimmed_seq(
            context_sequence_wt1_stranded,
            max(0, bp_in_wt1_nt_ex - self.context_seq_len),
            bp_in_wt1_nt_ex + self.context_seq_len
        )

        # the breakpoint in wt2 is at the position where the ft2 transcripts start
        bp_in_wt2_nt_ex = len(wt2_exon_transcripts) - len(ft2_exon_transcripts)
        context_sequence_wt2_stranded = get_stranded_seq(wt2_exon_transcripts, bp2.strand)
        context_sequence_wt2 = get_trimmed_seq(
            context_sequence_wt2_stranded,
            max(0, bp_in_wt2_nt_ex - self.context_seq_len),
            bp_in_wt2_nt_ex + self.context_seq_len
        )

        # Only contains the context sequence lengths not the position of the breakpoint
        context_sequence_bp = min(self.context_seq_len, bp_in_fusion_nt_ex)
        context_sequence_wt1_bp = min(self.context_seq_len, bp_in_wt1_nt_ex)
        context_sequence_wt2_bp = min(self.context_seq_len, bp_in_wt2_nt_ex)

        # == Protein AA sequence and Neo-peptide calculation ==

        # Create assumed fusion CDS sequence for peptide translation
        # TODO: What happens if there is no CDS in transcript 1?
        # Check if ft1's last CDS ends before the breakpoint. If yes, fusion will only contain the wildtype CDS.
        bp1_after_cds_end = False
        if bp1.strand == "+" and ft_cds_1 and ft_cds_1[-1].stop < bp1.pos:
            fusion_cds = ""
            bp1_after_cds_end = True
        elif bp1.strand == "-" and ft_cds_1 and ft_cds_1[0].start > bp1.pos:
            fusion_cds = ""
            bp1_after_cds_end = True
        # Classical in-frame fusion
        elif fusion_transcript.frame == "in_frame" and ft_cds_2:
            fusion_cds = get_stranded_seq(ft1_cds_transcripts, bp1.strand) + get_stranded_seq(ft2_cds_transcripts, bp2.strand)
        # In-frame fusion with no CDS on 2nd transcript (partial UTR fusion)
        elif fusion_transcript.frame == "in_frame" and not ft_cds_2:
            fusion_cds = get_stranded_seq(ft1_cds_transcripts, bp1.strand) + get_stranded_seq(ft2_exon_transcripts, bp2.strand)
            ft2_cds_transcripts = ft2_exon_transcripts
        # Not in-frame fusions, include full exon at fusion point from transcript 2
        elif fusion_transcript.frame != "in_frame":
            fusion_cds = get_stranded_seq(ft1_cds_transcripts, bp1.strand) + get_stranded_seq(ft2_exon_transcripts, bp2.strand)
            ft2_cds_transcripts = ft2_exon_transcripts
        else:
            fusion_cds = "NA"
        # Translate full fusion CDS
        # TODO: Use translation table 2 for mitochondrial genes
        translation_table = 1  # Standard table for human AA codes
        if fusion_cds in ["NA", ""]:
            fusion_peptide = ""
            fusion_protein_sequence = fusion_peptide
            fusion_protein_sequence_bp = 0
            neo_peptide_sequence = ""
            neo_peptide_sequence_bp = 0
        else:
            fusion_peptide = fusion_cds[wt1_full.frame_at_start:].translate(table=translation_table, to_stop=True)
            fusion_protein_sequence = fusion_peptide

            # Breakpoint in fusion transcript (CDS)
            bp_in_fusion_nt = len(ft1_cds_transcripts)

            # Breakpoint in fusion peptide, whereas x.0 / x.3 / x.6 indicate that the breakpoint 
            # is at the beginning / after first base / after second base of the underlying codon
            bp_in_fusion_aa = ((bp_in_fusion_nt - wt1_full.frame_at_start) * 10 // 3) / 10
            fusion_protein_sequence_bp = round(bp_in_fusion_aa, 1) # only used for printing, could be removed

            # Neo-peptide sequence breakpoint
            if (int(bp_in_fusion_aa) - 13) < 0:
                neo_peptide_sequence_bp = round(bp_in_fusion_aa, 1)
            else:
                neo_peptide_sequence_bp = round(bp_in_fusion_aa - (int(bp_in_fusion_aa) - 13), 1)

            # Neo-peptide
            neo_peptide_sequence = fusion_peptide[max(0, int(bp_in_fusion_aa) - 13) :]
            # Truncate the neo-peptide sequence to the breakpoint position if in-frame
            if fusion_transcript.frame == "in_frame":
                neo_peptide_sequence = neo_peptide_sequence[: (int(neo_peptide_sequence_bp) + 13)]

        # Check if this context sequence has already been created
        ftid = fusion_transcript.get_ftid()
        ftid_plus = f"{ftid}_{ctx_id}"
        if ftid_plus in self.ftid_plus_set:
            annotation_bias = True
        else:
            annotation_bias = False
            self.ftid_plus_set.add(ftid_plus)

        # Track edge cases for later on filtering. Default is "pass"
        filter_comments = []
        if bp1_after_cds_end:
            filter_comments.append("bp1 after cds end")
        # No CDS at the left fusion partner leads in both versions to no fusion peptide
        if not ft_cds_1:
            filter_comments.append("no ft1 cds")
        # No CDS at the right fusion partner in in-frame fusions may lead to incorrect fusion transcripts
        if not ft_cds_2 and fusion_transcript.frame == "in_frame":
            filter_comments.append("no ft2 cds")

        # Check if breakpoints are within exons of the fusion transcript
        # This could be optimized
        for exon in ft_exons_1:
            if exon.start <= bp1.pos <= exon.stop:
                break
        else:
            filter_comments.append("bp1 not in exons")
        for exon in ft_exons_2:
            if exon.start <= bp2.pos <= exon.stop:
                break
        else:
            filter_comments.append("bp2 not in exons")

        # Join all filter comments or use "pass" if none
        filter_comment = "|".join(filter_comments) if filter_comments else "pass"

        neo_pep_nuc_length = len(neo_peptide_sequence) * 3
        neo_pep_until_bp_nuc = int(round(neo_peptide_sequence_bp * 3, 0))

        # Where for is the exon starts and ends needed?
        # The result can include multiple exon starts and ends from both transcripts
        if bp1.strand == "+":
            exon_starts_1, exon_ends_1 = get_exon_ranges_reverse(ft_exons_1, neo_pep_until_bp_nuc)
        else:
            exon_starts_1, exon_ends_1 = get_exon_ranges(ft_exons_1, neo_pep_until_bp_nuc)

        if bp2.strand == "+":
            exon_starts_2, exon_ends_2 = get_exon_ranges(ft_exons_2, neo_pep_nuc_length - neo_pep_until_bp_nuc)
        else:
            exon_starts_2, exon_ends_2 = get_exon_ranges_reverse(ft_exons_2, neo_pep_nuc_length - neo_pep_until_bp_nuc)

        exon_starts = f"{exon_starts_1}*{exon_starts_2}"
        exon_ends = f"{exon_ends_1}*{exon_ends_2}"

        wt1_is_good_transcript = wt1_full.flags
        if len(wt1_cds_transcripts) % 3 != 0:
            wt1_is_good_transcript.add("wt1 seq % 3 != 0")
        wt2_is_good_transcript = wt2_full.flags
        if len(wt2_cds_transcripts) % 3 != 0:
            wt2_is_good_transcript.add("wt2 seq % 3 != 0")

        wt1_start_stop = ""
        wt2_start_stop = ""
        if wt1_full.exons:
            wt1_start_stop = f"{bp1.chrom}:{wt1_full.exons[0].start}:{wt1_full.exons[-1].stop}"
        if wt2_full.exons:
            wt2_start_stop = f"{bp2.chrom}:{wt2_full.exons[0].start}:{bp2.chrom}:{wt2_full.exons[-1].stop}"

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
            "ft1_exon_nr": str(len(ft_exons_1)),
            "ft2_exon_nr": str(len(wt2_full.exons) - len(ft_exons_2) + 1),
            "exon_starts": exon_starts,
            "exon_ends": exon_ends,
            "exon_boundary1": bp1.exon_boundary,
            "exon_boundary2": bp2.exon_boundary,
            "exon_boundary": fusion_transcript.get_combined_boundary(),
            "bp1_frame": str(wt1_full.frame_at_bp),
            "bp2_frame": str(wt2_full.frame_at_bp),
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
            "wt1_exon_pos": wt1_full.exons,
            "wt2_exon_pos": wt2_full.exons,
            "ft1_exon_pos": ft_exons_1,
            "ft2_exon_pos": ft_exons_2,
            "wt1_cds_pos": wt1_full.cds,
            "wt2_cds_pos": wt2_full.cds,
            "ft1_cds_pos": ft_cds_1,
            "ft2_cds_pos": ft_cds_2,
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
            "wt1_cds_transcripts": str(wt1_cds_transcripts),
            "wt2_cds_transcripts": str(wt2_cds_transcripts),
            "ft1_cds_transcripts": ft1_cds_transcripts,
            "ft2_cds_transcripts": ft2_cds_transcripts,
            "wt1_peptide": wt1_peptide,
            "wt2_peptide": wt2_peptide,
            "fusion_transcript": fusion_cds,  # should be renamed in the output
            "fusion_peptide": fusion_peptide,  # same as fusion_protein_sequence
            "wt1_is_good_transcript": wt1_is_good_transcript,
            "wt2_is_good_transcript": wt2_is_good_transcript,
            "wt1_trans_biotype": wt1_full.transcript_biotype,
            "wt2_trans_biotype": wt2_full.transcript_biotype,
            "wt1_gene_biotype": wt1_full.gene_biotype,
            "wt2_gene_biotype": wt2_full.gene_biotype,
            "wt1_description": wt1_full.description,
            "wt2_description": wt2_full.description,
            "wt1_frame_at_start": wt1_full.frame_at_start,
            "wt2_frame_at_start": wt2_full.frame_at_start,
            "wt1_TSL": wt1_full.tsl,
            "wt2_TSL": wt2_full.tsl,
            "wt1_exon_no": len(wt1_full.exons),
            "wt2_exon_no": len(wt2_full.exons),
            "ft1_exon_no": len(ft_exons_1),
            "ft2_exon_no": len(ft_exons_2),
            "wt1_cds_no": len(wt1_full.cds),
            "wt2_cds_no": len(wt2_full.cds),
            "ft1_cds_no": len(ft_cds_1),
            "ft2_cds_no": len(ft_cds_2),
            "wt1_start_stop": wt1_start_stop,
            "wt2_start_stop": wt2_start_stop,
            "annotation_bias": annotation_bias,
            "filter_comment": filter_comment
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
