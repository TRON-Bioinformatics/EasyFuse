"""
This module handles data loading.
"""

import csv

# pylint: disable=E0401
from breakpoint import Breakpoint

def load_detected_fusions(detected_fusions_table: str) -> dict:
    """Read BPIDs from Detected_Fusions.csv and store the info in a dict.

    Args:
        detected_fusions_table (str): Path to detected fusions table

    Returns:
        dict: Dictionary with BPIDs as keys and BPs as values, where
        BPIDs are the unique identifiers for each fusion event and BPs
        BPID: 21:41494380:-_7:13935843:-
        BP1: 21:41494380:-
        BP2: 7:13935843:-
    """

    bp_dict = {}
    with open(detected_fusions_table, "r", encoding="utf8") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        for row in reader:
            bpid = row["BPID"]
            bp1_str = row["Breakpoint1"]
            bp2_str = row["Breakpoint2"]
            bp1_chr, bp1_pos, bp1_strand = bp1_str.split(":")
            bp2_chr, bp2_pos, bp2_strand = bp2_str.split(":")
            bp1 = Breakpoint(bp1_chr, int(bp1_pos), bp1_strand)
            bp2 = Breakpoint(bp2_chr, int(bp2_pos), bp2_strand)
            # This could also be implemented as a set as BPID already contains both BPs
            if not bpid in bp_dict:
                bp_dict[bpid] = (bp1, bp2)
    return bp_dict


def load_tsl_data(tsl_in_file: str) -> dict:
    """Load data created by gtf2tsl.py into a dict"""
    #logger.info("Loading TSL data into dict")
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
