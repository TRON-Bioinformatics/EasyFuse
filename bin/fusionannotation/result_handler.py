

# pylint: disable=E0401
import xxhash # type: ignore
from Bio.Seq import Seq

from file_headers import FULL_ANNOTATION_HEADER
from sequence_validation import get_stranded_seq

header_idx = range(len(FULL_ANNOTATION_HEADER))
header_dict = dict(zip(FULL_ANNOTATION_HEADER, header_idx))

def update_results_list(results_lists, cds_seq_dict, context_seq_len):
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
    last_chr = ""
    ftid_plus_set = set()
    tmp_dict = None
    for result_list in results_lists:
        # create a tempory dict with the cds_pos tuple as keys and
        # the respective genomic sequence as value
        # the dict is restricted to the chromosome of the fusion partner and
        # only updated if the chromosome changes

        if not last_chr == result_list[header_dict["bp1_chr"]]:
            tmp_dict = dict(
                zip(
                    cds_seq_dict[result_list[header_dict["bp1_chr"]]][0],
                    cds_seq_dict[result_list[header_dict["bp1_chr"]]][1],
                )
            )
            last_chr = result_list[header_dict["bp1_chr"]]
        # add sequences as gffutils Seq objects to the table

        result_list[header_dict["wt1_exon_seqs"]] = []
        for feature in result_list[header_dict["wt1_exon_pos"]]:
            result_list[header_dict["wt1_exon_seqs"]].append(tmp_dict[feature])
        result_list[header_dict["ft1_exon_seqs"]] = [
            tmp_dict[feature]
            for feature in result_list[header_dict["ft1_exon_pos"]]
        ]
        result_list[header_dict["wt1_cds_seqs"]] = [
            tmp_dict[feature] for feature in result_list[header_dict["wt1_cds_pos"]]
        ]
        result_list[header_dict["ft1_cds_seqs"]] = [
            tmp_dict[feature] for feature in result_list[header_dict["ft1_cds_pos"]]
        ]
        # and concatenate them
        result_list[header_dict["wt1_exon_transcripts"]] = sum(
            result_list[header_dict["wt1_exon_seqs"]], Seq("")
        )
        result_list[header_dict["ft1_exon_transcripts"]] = sum(
            result_list[header_dict["ft1_exon_seqs"]], Seq("")
        )
        result_list[header_dict["wt1_cds_transcripts"]] = sum(
            result_list[header_dict["wt1_cds_seqs"]], Seq("")
        )
        result_list[header_dict["ft1_cds_transcripts"]] = sum(
            result_list[header_dict["ft1_cds_seqs"]], Seq("")
        )
        # create the first part of the context and fusion sequence and translate wt1 seqs
        translation_shift1 = result_list[
            header_dict["wt1_frame_at_start"]
        ]  # for easier visualization, store in separate var
        trans_table = 1
        if last_chr == "MT":
            trans_table = 2
        result_list[header_dict["context_sequence"]] = get_stranded_seq(
            result_list[header_dict["ft1_exon_transcripts"]],
            result_list[header_dict["bp1_strand"]],
        )
        result_list[header_dict["context_sequence_wt1"]] = get_stranded_seq(
            result_list[header_dict["wt1_exon_transcripts"]],
            result_list[header_dict["bp1_strand"]],
        )
        result_list[header_dict["fusion_transcript"]] = get_stranded_seq(
            result_list[header_dict["ft1_cds_transcripts"]],
            result_list[header_dict["bp1_strand"]],
        )
        result_list[header_dict["wt1_peptide"]] = get_stranded_seq(
            result_list[header_dict["wt1_cds_transcripts"]],
            result_list[header_dict["bp1_strand"]],
        )[translation_shift1:].translate(table=trans_table)

        # do the same for fusion partner 2
        if not last_chr == result_list[header_dict["bp2_chr"]]:
            tmp_dict = dict(
                zip(
                    cds_seq_dict[result_list[header_dict["bp2_chr"]]][0],
                    cds_seq_dict[result_list[header_dict["bp2_chr"]]][1],
                )
            )
            last_chr = result_list[header_dict["bp2_chr"]]
        result_list[header_dict["wt2_exon_seqs"]] = [
            tmp_dict[feature]
            for feature in result_list[header_dict["wt2_exon_pos"]]
        ]
        result_list[header_dict["ft2_exon_seqs"]] = [
            tmp_dict[feature]
            for feature in result_list[header_dict["ft2_exon_pos"]]
        ]
        result_list[header_dict["wt2_cds_seqs"]] = [
            tmp_dict[feature] for feature in result_list[header_dict["wt2_cds_pos"]]
        ]
        result_list[header_dict["ft2_cds_seqs"]] = [
            tmp_dict[feature] for feature in result_list[header_dict["ft2_cds_pos"]]
        ]
        result_list[header_dict["wt2_exon_transcripts"]] = sum(
            result_list[header_dict["wt2_exon_seqs"]], Seq("")
        )
        result_list[header_dict["ft2_exon_transcripts"]] = sum(
            result_list[header_dict["ft2_exon_seqs"]], Seq("")
        )
        result_list[header_dict["wt2_cds_transcripts"]] = sum(
            result_list[header_dict["wt2_cds_seqs"]], Seq("")
        )
        result_list[header_dict["ft2_cds_transcripts"]] = sum(
            result_list[header_dict["ft2_cds_seqs"]], Seq("")
        )
        translation_shift2 = result_list[header_dict["wt2_frame_at_start"]]
        # append the second part of the context and fusion sequence and translate wt2 seqs
        trans_table = 1
        if last_chr == "MT":
            trans_table = 2
        result_list[header_dict["context_sequence"]] += get_stranded_seq(
            result_list[header_dict["ft2_exon_transcripts"]],
            result_list[header_dict["bp2_strand"]],
        )
        result_list[header_dict["context_sequence_wt2"]] = get_stranded_seq(
            result_list[header_dict["wt2_exon_transcripts"]],
            result_list[header_dict["bp2_strand"]],
        )
        # for in frame fusions, use the cds sequences, for everything else, use the exons
        if result_list[header_dict["frame"]] == "in_frame":
            result_list[header_dict["fusion_transcript"]] += get_stranded_seq(
                result_list[header_dict["ft2_cds_transcripts"]],
                result_list[header_dict["bp2_strand"]],
            )
        else:
            result_list[header_dict["fusion_transcript"]] += get_stranded_seq(
                result_list[header_dict["ft2_exon_transcripts"]],
                result_list[header_dict["bp2_strand"]],
            )
        result_list[header_dict["wt2_peptide"]] = get_stranded_seq(
            result_list[header_dict["wt2_cds_transcripts"]],
            result_list[header_dict["bp2_strand"]],
        )[translation_shift2:].translate(table=trans_table)

        # bp pos in the fusion
        # the fusion starts in the fusion transcript with the end of the first part
        bp_in_fusion_nt = len(result_list[header_dict["ft1_cds_transcripts"]])
        # the fusion starts in the fusion transcript with the end of the first part
        bp_in_fusion_nt_ex = len(result_list[header_dict["ft1_exon_transcripts"]])
        # the breakpoint in wt1 must be identical to the breakpoint in the fusion
        bp_in_wt1_nt_ex = bp_in_fusion_nt_ex
        # the breakpoint in wt2 is at the position where the ft2 transcripts start
        bp_in_wt2_nt_ex = len(
            result_list[header_dict["wt2_exon_transcripts"]]
        ) - len(result_list[header_dict["ft2_exon_transcripts"]])
        # the respective bp pos in the aa seq, whereas x.0 / x.3 / x.6
        # indicate that the breakpoint is at the beginning / after first base /
        # after second base of the underlying codon
        bp_in_fusion_aa = (
            (bp_in_fusion_nt - translation_shift1) * 10.0 // 3
        ) / 10.0

        # translate and trim the fusion sequence around the breakpoint
        result_list[header_dict["fusion_peptide"]] = result_list[
            header_dict["fusion_transcript"]
        ][translation_shift1:].translate(table=1, to_stop=True)
        result_list[header_dict["context_sequence_100"]] = result_list[
            header_dict["context_sequence"]
        ][max(0, bp_in_fusion_nt_ex - 100) : (bp_in_fusion_nt_ex + 100)]
        result_list[header_dict["context_sequence"]] = result_list[
            header_dict["context_sequence"]
        ][
            max(0, bp_in_fusion_nt_ex - context_seq_len) : (
                bp_in_fusion_nt_ex + context_seq_len
            )
        ]
        result_list[header_dict["context_sequence_wt1"]] = result_list[
            header_dict["context_sequence_wt1"]
        ][
            max(0, bp_in_wt1_nt_ex - context_seq_len) : (
                bp_in_wt1_nt_ex + context_seq_len
            )
        ]
        result_list[header_dict["context_sequence_wt2"]] = result_list[
            header_dict["context_sequence_wt2"]
        ][
            max(0, bp_in_wt2_nt_ex - context_seq_len) : (
                bp_in_wt2_nt_ex + context_seq_len
            )
        ]
        result_list[header_dict["context_sequence_bp"]] = min(
            context_seq_len, bp_in_fusion_nt_ex
        )
        result_list[header_dict["context_sequence_wt1_bp"]] = min(
            context_seq_len, bp_in_wt1_nt_ex
        )
        result_list[header_dict["context_sequence_wt2_bp"]] = min(
            context_seq_len, bp_in_wt2_nt_ex
        )

        result_list[header_dict["neo_peptide_sequence"]] = result_list[
            header_dict["fusion_peptide"]
        ][max(0, int(bp_in_fusion_aa) - 13) :]

        result_list[header_dict["fusion_protein_sequence"]] = result_list[
            header_dict["fusion_peptide"]
        ]

        result_list[header_dict["fusion_protein_sequence_bp"]] = round(
                bp_in_fusion_aa, 1
            )

        if (int(bp_in_fusion_aa) - 13) < 0:
            result_list[header_dict["neo_peptide_sequence_bp"]] = round(
                bp_in_fusion_aa, 1
            )
        else:
            result_list[header_dict["neo_peptide_sequence_bp"]] = round(
                bp_in_fusion_aa - (int(bp_in_fusion_aa) - 13), 1
            )
        if result_list[header_dict["frame"]] == "in_frame":
            result_list[header_dict["neo_peptide_sequence"]] = result_list[
                header_dict["neo_peptide_sequence"]
            ][: (int(result_list[header_dict["neo_peptide_sequence_bp"]]) + 13)]

        result_list[header_dict["context_sequence_id"]] = xxhash.xxh64(
            str(result_list[header_dict["context_sequence"]])
        ).hexdigest()
        result_list[header_dict["context_sequence_100_id"]] = xxhash.xxh64(
            str(result_list[header_dict["context_sequence_100"]])
        ).hexdigest()

        ftid = result_list[header_dict["FTID"]]
        ctx_id = result_list[header_dict["context_sequence_id"]]
        ftid_plus = f"{ftid}_{ctx_id}"
        if ftid_plus in ftid_plus_set:
            result_list[header_dict["annotation_bias"]] = True
        else:
            result_list[header_dict["annotation_bias"]] = False
            ftid_plus_set.add(ftid_plus)

        # set exon starts/ends
        # urla: needs revision and simplification, but is working...
        neo_pep_nuc_length = (
            len(result_list[header_dict["neo_peptide_sequence"]]) * 3
        )
        neo_pep_until_bp_nuc = int(
            round(result_list[header_dict["neo_peptide_sequence_bp"]] * 3, 0)
        )
        exon_pos = 1
        if result_list[header_dict["bp1_strand"]] == "+":
            summed_len = 0
            for exon_pos, exon_length in enumerate(
                [y - x for x, y in result_list[header_dict["ft1_exon_pos"]]][::-1],
                1,
            ):
                summed_len += exon_length
                if neo_pep_until_bp_nuc <= summed_len:
                    break
            lookup_exons1 = result_list[header_dict["ft1_exon_pos"]][-exon_pos:]
        else:
            summed_len = 0
            for exon_pos, exon_length in enumerate(
                [y - x for x, y in result_list[header_dict["ft1_exon_pos"]]], 1
            ):
                summed_len += exon_length
                if neo_pep_until_bp_nuc <= summed_len:
                    break
            lookup_exons1 = result_list[header_dict["ft1_exon_pos"]][:exon_pos]

        if result_list[header_dict["bp2_strand"]] == "+":
            summed_len = 0
            for exon_pos, exon_length in enumerate(
                [y - x for x, y in result_list[header_dict["ft2_exon_pos"]]], 1
            ):
                summed_len += exon_length
                if (neo_pep_nuc_length - neo_pep_until_bp_nuc) <= summed_len:
                    break
            lookup_exons2 = result_list[header_dict["ft2_exon_pos"]][:exon_pos]
        else:
            summed_len = 0
            for exon_pos, exon_length in enumerate(
                [y - x for x, y in result_list[header_dict["ft2_exon_pos"]]][::-1],
                1,
            ):
                summed_len += exon_length
                if (neo_pep_nuc_length - neo_pep_until_bp_nuc) <= summed_len:
                    break
            lookup_exons2 = result_list[header_dict["ft2_exon_pos"]][-exon_pos:]

        exon_starts1_str = "|".join([str(x) for x, y in lookup_exons1])
        exon_starts2_str = "|".join([str(x) for x, y in lookup_exons2])
        result_list[header_dict["exon_starts"]] = f"{exon_starts1_str}*{exon_starts2_str}"

        exon_ends1_str = "|".join([str(y) for x, y in lookup_exons1])
        exon_ends2_str = "|".join([str(y) for x, y in lookup_exons2])
        result_list[header_dict["exon_ends"]] = f"{exon_ends1_str}*{exon_ends2_str}"

        # experimental
        if not len(result_list[header_dict["wt1_cds_transcripts"]]) % 3 == 0:
            result_list[header_dict["wt1_is_good_transcript"]].add(
                "wt1 seq % 3 != 0"
            )
        if not len(result_list[header_dict["wt2_cds_transcripts"]]) % 3 == 0:
            result_list[header_dict["wt2_is_good_transcript"]].add(
                "wt2 seq % 3 != 0"
            )

        result_list[header_dict["wt1_cds_transcripts"]] = str(
            result_list[header_dict["wt1_cds_transcripts"]]
        )
        result_list[header_dict["wt2_cds_transcripts"]] = str(
            result_list[header_dict["wt2_cds_transcripts"]]
        )
    return results_lists
