"""
This module provides methods for sequence operations.
"""

# pylint: disable=import-error
import xxhash # type: ignore
from Bio.Seq import Seq # type: ignore

from bin.fusionannotation.src.breakpoint import Breakpoint
from bin.fusionannotation.src.fusion_transcript import FusionTranscript


def get_stranded_seq(sequence: Seq, strand: str) -> Seq:
    """
    Return the reverse complement of a sequence 
    if the strand is negative and the unchanged sequence otherwise.
    """
    if strand == "-":
        return sequence.reverse_complement()
    return sequence


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
        ft1_seq (Seq): Sequence of the first fusion transcript
        f2_seq (Seq): Sequence of the second fusion transcript
        bp1_strand (str): Strand of the first breakpoint
        bp2_strand (str): Strand of the second breakpoint

    Returns:
        Seq: Context sequence
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
    """Returns the fusion transcript sequence.

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


def get_peptide_sequence(transcript_sequence: Seq, bp: Breakpoint, transcript_frame: int) -> Seq:
    """Returns the peptide sequence for a transcript.

    Args:
        transcript_sequence (Seq): Sequence to be translated
        bp (Breakpoint): Breakpoint object
        transcript_frame (int): Frame shift of the transcript

    Returns:
        Seq: Peptide sequence
    """
    # for easier visualization, store in separate var
    trans_table = 2 if bp.chrom == "MT" else 1

    peptide_seq_raw = get_stranded_seq(
        transcript_sequence,
        bp.strand
    )
    peptide_seq = peptide_seq_raw[transcript_frame:].translate(table=trans_table)
    return peptide_seq


def calc_hash(seq: Seq) -> str:
    """Calculates the hash of a sequence.

    Args:
        seq (Seq): Sequence to calculate the hash for

    Returns:
        str: Hash of the sequence
    """
    return xxhash.xxh64(str(seq)).hexdigest()
