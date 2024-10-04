"""
This module handles data loading.
"""


def breakpoints_to_dict(bp_in_file):
    """Read BPIDs from Detected_Fusions.csv and store the info in a dict"""
    #logger.info("Loading BPIDs from Detected_Fusions.csv")
    bp_dict = {}
    with open(bp_in_file, "r", encoding="utf8") as bpin:
        next(bpin)  # skip header
        for line in bpin:
            elements = line.rstrip().split(";")
            bpid = elements[0]
            bp1 = elements[2]
            bp2 = elements[3]
            # This could also be implemented as a set as BPID already contains both BPs
            if not bpid in bp_dict:
                bp_dict[bpid] = (bp1, bp2)
    return bp_dict


def load_tsl_data(tsl_in_file):
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
