"""
This module handles data loading.
"""

import csv

# pylint: disable=E0401
from Bio import SeqIO # type: ignore

from .breakpoint import Breakpoint

def load_detected_fusions(detected_fusions_table: str) -> dict:
    """Read BPIDs from Detected_Fusions.csv and store the info in a dict.

    Args:
        detected_fusions_table (str): Path to detected fusions table

    Returns:
        dict: Dictionary with BPIDs (e.g. 21:41494380:-_7:13935843:-) as keys
        and breakpoints (BP1, BP2) as values
    """

    bp_dict = {}
    with open(detected_fusions_table, "r", encoding="utf8") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            bp1_chr, bp1_pos, bp1_strand = row["Breakpoint1"].split(":")
            bp2_chr, bp2_pos, bp2_strand = row["Breakpoint2"].split(":")
            bp1 = Breakpoint(bp1_chr, int(bp1_pos), bp1_strand)
            bp2 = Breakpoint(bp2_chr, int(bp2_pos), bp2_strand)
            # This could also be implemented as a set as BPID already contains both BPs
            if not row["BPID"] in bp_dict:
                bp_dict[row["BPID"]] = (bp1, bp2)
    return bp_dict


def load_tsl_data(tsl_in_file: str) -> dict:
    """Load data created by gtf2tsl.py into a dict"""
    tsl_dict = {}
    with open(tsl_in_file, "r", encoding="utf8") as tslin:
        next(tslin)  # skip header
        for line in tslin:
            elements = line.rstrip().split("\t")
            trans_id = elements[0]
            tsl = elements[3]
            if not trans_id in tsl_dict:
                tsl_dict[trans_id] = tsl
    return tsl_dict


def load_genomic_data(genome_fasta: str, fusion_transcripts: list) -> dict:
    """Load genomic data from genome fasta and fusion transcripts"""
    feature_seqs = {}
    for ft in fusion_transcripts:
        chr1 = ft.bp1.chrom
        chr2 = ft.bp2.chrom
        for chrom in (chr1, chr2):
            if not chrom in feature_seqs:
                feature_seqs[chrom] = [set(), set()]
        for features in (
            ft.transcript_1.exons,
            ft.exons_transcript_1,
            ft.transcript_1.cds,
            ft.cds_transcript_1
        ):
            feature_seqs[chr1][0].update(features)
        for features in (
            ft.transcript_2.exons,
            ft.exons_transcript_2,
            ft.transcript_2.cds,
            ft.cds_transcript_2
        ):
            feature_seqs[chr2][0].update(features)
    for record in SeqIO.parse(genome_fasta, "fasta"):
        if record.id in feature_seqs:
            feature_seqs[record.id][1] = [
                record.seq[max(0, feature.start - 1) : feature.stop]
                for feature in feature_seqs[record.id][0]
            ]
    return feature_seqs
