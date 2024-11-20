

# pylint: disable=E0401
import xxhash # type: ignore
from Bio.Seq import Seq # type: ignore

from bin.fusionannotation.src.breakpoint import Breakpoint
from bin.fusionannotation.src.fusion_transcript import FusionTranscript
from bin.fusionannotation.src.sequence_validation import get_stranded_seq


def get_feature_seqs(tmp_dict: dict, features: list) -> list:
    """
    Return a list of Seq objects for a list of features.
    """
    return [tmp_dict[feature] for feature in features]


def concatenate_seqs(sequences: list) -> str:
    """Returns a concatenated string of Seq objects.
    If sequences is empty, an empty Seq is returned.

    Args:
        sequences (list): List of Seq objects

    Returns:
        Seq: Concatenated Seq object, empty if sequences is empty
    """
    return sum(sequences, Seq(""))


def get_context_sequence(ft1_seq: Seq, ft2_seq: Seq, bp1_strand: str, bp2_strand: str) -> str:
    """Returns the context sequence for a fusion transcript.

    Args:
        transcripts (Seq): _description_
        fusion_transcript (FusionTranscript): _description_

    Returns:
        str: _description_
    """
    ctx_part_1 = get_stranded_seq(ft1_seq, bp1_strand)
    ctx_part_2 = get_stranded_seq(ft2_seq, bp2_strand)
    return ctx_part_1 + ctx_part_2


def get_fusion_transcript_sequence(
    bp1: Breakpoint,
    bp2: Breakpoint,
    table_row: dict,
    fusion_transcript: FusionTranscript
) -> str:
    """_summary_

    Args:
        bp1 (Breakpoint): _description_
        bp2 (Breakpoint): _description_
        table_row (dict): _description_
        fusion_transcript (FusionTranscript): _description_

    Returns:
        str: _description_
    """
    fusion_transcript_sequence = get_stranded_seq(
        table_row["ft1_cds_transcripts"],
        bp1.strand
    )

    # for in frame fusions, use the cds sequences, for everything else, use the exons
    if fusion_transcript.frame == "in_frame":
        fusion_transcript_sequence += get_stranded_seq(
            fusion_transcript.cds_transcript_2,
            bp2.strand
        )
    else:
        fusion_transcript_sequence += get_stranded_seq(
            table_row["ft2_exon_transcripts"],
            bp2.strand
        )
    return fusion_transcript_sequence


def get_peptide_sequence(cds_transcripts: str, bp: object, transcript: object) -> str:
    """create the first part of the context and fusion sequence and translate wt1 seqs"""
    translation_shift = transcript.frame
    # for easier visualization, store in separate var
    trans_table = 2 if bp.chrom == "MT" else 1

    peptide_seq_raw = get_stranded_seq(
        cds_transcripts,
        bp.strand
    )
    peptide_seq = peptide_seq_raw[translation_shift:].translate(table=trans_table)
    return peptide_seq


def calc_hash(seq: Seq) -> str:
    """_summary_

    Args:
        seq (Seq): _description_

    Returns:
        str: _description_
    """
    return xxhash.xxh64(str(seq)).hexdigest()


def generate_temp_dict(cds_seq_dict: dict, chrom: str) -> dict:
    """_summary_

    Args:
        cds_seq_dict (dict): _description_
        bp (object): _description_

    Returns:
        dict: _description_
    """
    return dict(
        zip(
            cds_seq_dict[chrom][0],
            cds_seq_dict[chrom][1],
        )
    )


def get_exon_ranges(exons: list, neo_pep_until_bp_nuc: int) -> list:
    """_summary_

    Args:
        exons (list): _description_
        neo_pep_until_bp_nuc (int): _description_

    Returns:
        list: _description_
    """
    exon_pos = 1
    summed_len = 0
    for exon_pos, exon_length in enumerate([y - x for x, y in exons], 1):
        summed_len += exon_length
        if neo_pep_until_bp_nuc <= summed_len:
            break
    lookup_exons = exons[:exon_pos]
    exon_starts = "|".join([str(x) for x, y in lookup_exons])
    exon_ends = "|".join([str(y) for x, y in lookup_exons])
    return (exon_starts, exon_ends)


def get_exon_ranges_reverse(exons: list, neo_pep_until_bp_nuc: int) -> tuple:
    """_summary_

    Args:
        exons (list): _description_
        neo_pep_until_bp_nuc (int): _description_

    Returns:
        tuple: _description_
    """
    exon_pos = 1
    summed_len = 0
    for exon_pos, exon_length in enumerate([y - x for x, y in exons][::-1], 1):
        summed_len += exon_length
        if neo_pep_until_bp_nuc <= summed_len:
            break
    lookup_exons = exons[-exon_pos:]
    exon_starts = "|".join([str(x) for x, y in lookup_exons])
    exon_ends = "|".join([str(y) for x, y in lookup_exons])
    return (exon_starts, exon_ends)


def update_results_list(fusion_transcripts: list, cds_seq_dict: dict, context_seq_len: int) -> list:
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
    ftid_plus_set = set()
    result_list = []
    for fusion_transcript in fusion_transcripts:
        # create a tempory dict with the cds_pos tuple as keys and
        # the respective genomic sequence as value
        # the dict is restricted to the chromosome of the fusion partner and
        # only updated if the chromosome changes

        wt1 = fusion_transcript.transcript_1
        wt2 = fusion_transcript.transcript_2
        bp1 = fusion_transcript.bp1
        bp2 = fusion_transcript.bp2
        tmp_dict_1 = generate_temp_dict(cds_seq_dict, bp1.chrom)
        tmp_dict_2 = generate_temp_dict(cds_seq_dict, bp2.chrom)
        table_row = {
            "BPID": fusion_transcript.get_bpid(),
            "Fusion_Gene": fusion_transcript.get_fusion_gene_name(),
            "Breakpoint1": bp1,
            "Breakpoint2": bp2,
            "FTID": fusion_transcript.get_ftid(),
            "context_sequence_id": "context_sequence_id",
            "context_sequence_100_id": "context_sequence_100_id",
            "type": fusion_transcript.fusion_type,
            "exon_nr": str(fusion_transcript.get_exon_nr()),
            "ft1_exon_nr": str(len(fusion_transcript.exons_transcript_1)),
            "ft2_exon_nr": str(len(wt1.exon_pos_list) - len(fusion_transcript.exons_transcript_2)),
            "exon_starts": "exon_starts",
            "exon_ends": "exon_ends",
            "exon_boundary1": bp1.exon_boundary,
            "exon_boundary2": bp2.exon_boundary,
            "exon_boundary": fusion_transcript.get_combined_boundary(),
            "bp1_frame": str(wt1.frame),
            "bp2_frame": str(wt2.frame),
            "frame": fusion_transcript.frame,
            "context_sequence": "context_sequence",
            "context_sequence_bp": "context_sequence_bp",
            "neo_peptide_sequence": "neo_peptide_sequence",
            "neo_peptide_sequence_bp": "neo_peptide_sequence_bp",
            "fusion_protein_sequence": "fusion_protein_sequence",
            "fusion_protein_sequence_bp": "fusion_protein_sequence_bp",
            "context_sequence_wt1": "context_sequence_wt1",
            "context_sequence_wt2": "context_sequence_wt2",
            "context_sequence_wt1_bp": "context_sequence_wt1_bp",
            "context_sequence_wt2_bp": "context_sequence_wt2_bp",
            "context_sequence_100": "context_sequence_100",
            "bp1_chr": bp1.chrom,
            "bp1_pos": bp1.pos,
            "bp1_strand": bp1.strand,
            "bp2_chr": bp2.chrom,
            "bp2_pos": bp2.pos,
            "bp2_strand": bp2.strand,
            "cds_boundary1": bp1.cds_boundary,
            "cds_boundary2": bp2.cds_boundary,
            "wt1_exon_pos": wt1.exon_pos_list,
            "wt2_exon_pos": wt2.exon_pos_list,
            "ft1_exon_pos": fusion_transcript.exons_transcript_1,
            "ft2_exon_pos": fusion_transcript.exons_transcript_2,
            "wt1_cds_pos": wt1.cds_pos_list,
            "wt2_cds_pos": wt2.cds_pos_list,
            "ft1_cds_pos": fusion_transcript.cds_transcript_1,
            "ft2_cds_pos": fusion_transcript.cds_transcript_2,
            "wt1_exon_seqs": [],
            "wt2_exon_seqs": [],
            "ft1_exon_seqs": [],
            "ft2_exon_seqs": [],
            "wt1_cds_seqs": [],
            "wt2_cds_seqs": [],
            "ft1_cds_seqs": [],
            "ft2_cds_seqs": [],
            "wt1_exon_transcripts": "wt1_exon_transcripts",
            "wt2_exon_transcripts": "wt2_exon_transcripts",
            "ft1_exon_transcripts": "ft1_exon_transcripts",
            "ft2_exon_transcripts": "ft2_exon_transcripts",
            "wt1_cds_transcripts": "wt1_cds_transcripts",
            "wt2_cds_transcripts": "wt2_cds_transcripts",
            "ft1_cds_transcripts": "ft1_cds_transcripts",
            "ft2_cds_transcripts": "ft2_cds_transcripts",
            "wt1_peptide": "wt1_peptide",
            "wt2_peptide": "wt2_peptide",
            "fusion_transcript": "fusion_transcript",
            "fusion_peptide": "fusion_peptide",
            "wt1_is_good_transcript": wt1.is_good_transcript,
            "wt2_is_good_transcript": wt2.is_good_transcript,
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
            "wt1_exon_no": len(wt1.exon_pos_list),
            "wt2_exon_no": len(wt2.exon_pos_list),
            "ft1_exon_no": len(fusion_transcript.exons_transcript_1),
            "ft2_exon_no": len(fusion_transcript.exons_transcript_2),
            "wt1_cds_no": len(wt1.cds_pos_list),
            "wt2_cds_no": len(wt2.cds_pos_list),
            "ft1_cds_no": len(fusion_transcript.cds_transcript_1),
            "ft2_cds_no": len(fusion_transcript.cds_transcript_2),
            "wt1_start_stop": f"{bp1.chrom}:{wt1.exon_pos_list[0][0]}:{wt1.exon_pos_list[-1][1]}",
            "wt2_start_stop": f"{bp2.chrom}:{wt2.exon_pos_list[0][0]}:{wt2.exon_pos_list[-1][1]}",
            "annotation_bias": "annotation_bias",
        }

        table_row["wt1_cds_seqs"] = get_feature_seqs(tmp_dict_1, wt1.cds_pos_list)
        table_row["wt2_cds_seqs"] = get_feature_seqs(tmp_dict_2, wt2.cds_pos_list)
        table_row["wt1_cds_transcripts"] = concatenate_seqs(table_row["wt1_cds_seqs"])
        table_row["wt2_cds_transcripts"] = concatenate_seqs(table_row["wt2_cds_seqs"])

        table_row["ft1_cds_seqs"] = get_feature_seqs(tmp_dict_1, fusion_transcript.cds_transcript_1)
        table_row["ft2_cds_seqs"] = get_feature_seqs(tmp_dict_2, fusion_transcript.cds_transcript_2)
        table_row["ft1_cds_transcripts"] = concatenate_seqs(table_row["ft1_cds_seqs"])
        table_row["ft2_cds_transcripts"] = concatenate_seqs(table_row["ft2_cds_seqs"])

        table_row["ft1_exon_seqs"] = get_feature_seqs(tmp_dict_1, fusion_transcript.exons_transcript_1)
        table_row["ft2_exon_seqs"] = get_feature_seqs(tmp_dict_2, fusion_transcript.exons_transcript_2)
        table_row["ft1_exon_transcripts"] = concatenate_seqs(table_row["ft1_exon_seqs"])
        table_row["ft2_exon_transcripts"] = concatenate_seqs(table_row["ft2_exon_seqs"])

        table_row["wt1_exon_seqs"] = get_feature_seqs(tmp_dict_1, wt1.exon_pos_list)
        table_row["wt2_exon_seqs"] = get_feature_seqs(tmp_dict_2, wt2.exon_pos_list)
        table_row["wt1_exon_transcripts"] = concatenate_seqs(table_row["wt1_exon_seqs"])
        table_row["wt2_exon_transcripts"] = concatenate_seqs(table_row["wt2_exon_seqs"])

        table_row["wt1_peptide"] = get_peptide_sequence(table_row["wt1_cds_transcripts"], bp1, wt1)
        table_row["wt2_peptide"] = get_peptide_sequence(table_row["wt2_cds_transcripts"], bp2, wt2)

        table_row["fusion_transcript"] = get_fusion_transcript_sequence(
            bp1, bp2, table_row, fusion_transcript
        )

        # translate and trim the fusion sequence around the breakpoint
        table_row["fusion_peptide"] = table_row["fusion_transcript"][fusion_transcript.frame:].translate(table=1, to_stop=True)

        # the fusion starts in the fusion transcript with the end of the first part
        bp_in_fusion_nt_ex = len(fusion_transcript.exon_transcript_2)
        table_row["context_sequence_100"] = table_row["context_sequence"][
            max(0, bp_in_fusion_nt_ex - 100) : (bp_in_fusion_nt_ex + 100)
        ]
        table_row["context_sequence"] = get_context_sequence(
            table_row["ft1_exon_transcripts"],
            table_row["ft2_exon_transcripts"],
            bp1.strand,
            bp2.strand
        )[
            max(0, bp_in_fusion_nt_ex - context_seq_len) : (
                bp_in_fusion_nt_ex + context_seq_len
            )
        ]

        table_row["context_sequence_id"] = calc_hash(table_row["context_sequence"])
        table_row["context_sequence_100_id"] = calc_hash(table_row["context_sequence_100"])

        # the breakpoint in wt1 must be identical to the breakpoint in the fusion
        bp_in_wt1_nt_ex = bp_in_fusion_nt_ex
        table_row["context_sequence_wt1"] = get_stranded_seq(
            table_row["wt1_exon_transcripts"],
            bp1.strand
        )[
            max(0, bp_in_wt1_nt_ex - context_seq_len) : (
                bp_in_wt1_nt_ex + context_seq_len
            )
        ]

        # the breakpoint in wt2 is at the position where the ft2 transcripts start
        bp_in_wt2_nt_ex = len(table_row["wt2_exon_transcripts"]) - len(table_row["ft2_exon_transcripts"])
        table_row["context_sequence_wt2"] = get_stranded_seq(
            wt2.exon_pos_list,
            bp2.strand
        )[
            max(0, bp_in_wt2_nt_ex - context_seq_len) : (
                bp_in_wt2_nt_ex + context_seq_len
            )
        ]
        table_row["context_sequence_bp"] = min(context_seq_len, bp_in_fusion_nt_ex)
        table_row["context_sequence_wt1_bp"] = min(context_seq_len, bp_in_wt1_nt_ex)
        table_row["context_sequence_wt2_bp"] = min(context_seq_len, bp_in_wt2_nt_ex)

        # bp pos in the fusion
        # the fusion starts in the fusion transcript with the end of the first part
        bp_in_fusion_nt = len(fusion_transcript.cds_transcript_2)
        # the respective bp pos in the aa seq, whereas x.0 / x.3 / x.6
        # indicate that the breakpoint is at the beginning / after first base /
        # after second base of the underlying codon
        bp_in_fusion_aa = (
            (bp_in_fusion_nt - fusion_transcript.frame) * 10.0 // 3
        ) / 10.0
        table_row["neo_peptide_sequence"] = table_row["fusion_peptide"][
            max(0, int(bp_in_fusion_aa) - 13) :
        ]

        table_row["fusion_protein_sequence"] = table_row["fusion_peptide"]

        table_row["fusion_protein_sequence_bp"] = round(bp_in_fusion_aa, 1)

        if (int(bp_in_fusion_aa) - 13) < 0:
            table_row["neo_peptide_sequence_bp"] = round(bp_in_fusion_aa, 1)
        else:
            table_row["neo_peptide_sequence_bp"] = round(
                bp_in_fusion_aa - (int(bp_in_fusion_aa) - 13), 1
            )
        if fusion_transcript.frame == "in_frame":
            table_row["neo_peptide_sequence"] = table_row["neo_peptide_sequence"][
                : (int(table_row["neo_peptide_sequence_bp"]) + 13)
            ]


        ftid = table_row["FTID"]
        ctx_id = table_row["context_sequence_id"]
        ftid_plus = f"{ftid}_{ctx_id}"
        if ftid_plus in ftid_plus_set:
            table_row["annotation_bias"] = True
        else:
            table_row["annotation_bias"] = False
            ftid_plus_set.add(ftid_plus)

        # set exon starts/ends
        # urla: needs revision and simplification, but is working...
        neo_pep_nuc_length = (
            len(table_row["neo_peptide_sequence"]) * 3
        )
        neo_pep_until_bp_nuc = int(
            round(table_row["neo_peptide_sequence_bp"] * 3, 0)
        )

        if bp1.strand == "+":
            (exon_starts_1, exon_ends_1) = get_exon_ranges_reverse(
                table_row["ft1_exon_pos"],
                neo_pep_until_bp_nuc
            )
        else:
            (exon_starts_1, exon_ends_1) = get_exon_ranges(
                table_row["ft1_exon_pos"],
                neo_pep_until_bp_nuc
            )

        if bp2.strand == "+":
            (exon_starts_2, exon_ends_2) = get_exon_ranges(
                table_row["ft2_exon_pos"],
                neo_pep_nuc_length - neo_pep_until_bp_nuc
            )
        else:
            (exon_starts_2, exon_ends_2) = get_exon_ranges_reverse(
                table_row["ft2_exon_pos"],
                neo_pep_nuc_length - neo_pep_until_bp_nuc
            )

        table_row["exon_starts"] = f"{exon_starts_1}*{exon_starts_2}"
        table_row["exon_ends"] = f"{exon_ends_1}*{exon_ends_2}"

        # experimental
        if len(table_row["wt1_cds_transcripts"]) % 3 != 0:
            table_row["wt1_is_good_transcript"].add(
                "wt1 seq % 3 != 0"
            )
        if len(table_row["wt2_cds_transcripts"]) % 3 != 0:
            table_row["wt2_is_good_transcript"].add(
                "wt2 seq % 3 != 0"
            )

        table_row["wt1_cds_transcripts"] = str(
            table_row["wt1_cds_transcripts"]
        )
        table_row["wt2_cds_transcripts"] = str(
            table_row["wt2_cds_transcripts"]
        )
    return result_list
