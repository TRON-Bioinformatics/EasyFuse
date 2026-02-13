"""
Module for parsing infusion results
"""

import csv
import operator

from .file_headers import CHROMOSOMES


def get_fusion_gene_id_infusion(gene_id: str, transcript_list_field: str) -> str:
    """
    Helper method for infusion data parsing.
    Returns the most frequently listed fusion
    id in a list of fusion ids
    """
    # if only 1 gene id is returned as fusion partner, return this
    if not ";" in gene_id:
        return gene_id.replace("_", "-")
    # if 2 or more gene ids are given, take the most frequent from the transcript list
    gid_dict = {}
    for gid in transcript_list_field.split(";"):
        id_test = gid.split(":")[0][:-4]
        if not id_test in gid_dict:
            gid_dict[id_test] = 0
        gid_dict[id_test] += 1

    best_hit = max(iter(gid_dict.items()), key=operator.itemgetter(1))[0]

    return best_hit.replace("_", "-")


def parse_infusion_results(infile: str) -> list:
    """
    Load and parse results from infusion
    Description of output here: https://bitbucket.org/kokonech/infusion/wiki/Home
    """
    fusions = []
    with open(infile, encoding="utf8") as prediction:
        reader = csv.DictReader(prediction, delimiter="\t")
        for row in reader:
            up_gene = get_fusion_gene_id_infusion(row["gene_1"], row["transcript_1"])
            dn_gene = get_fusion_gene_id_infusion(row["gene_2"], row["transcript_2"])
            fusion_gene = f"{up_gene}_{dn_gene}".upper()

            # check whether fusion gene is not on primary chr
            if row["ref1"] not in CHROMOSOMES or row["ref2"] not in CHROMOSOMES:
                continue

            # element[1/4] = chr num, [2/5] = breakpoints, [23/29] = strands
            up_gene_bp = row["ref1"] + ":" + row["break_pos1"] + ":" + row["gene_1_strand"]
            dn_gene_bp = row["ref2"] + ":" + row["break_pos2"] + ":" + row["gene_2_strand"]

            bpid = up_gene_bp + "_" + dn_gene_bp

            fusions.append(
                {
                    "BPID": bpid,
                    "Fusion_Gene": fusion_gene,
                    "Breakpoint1": up_gene_bp,
                    "Breakpoint2": dn_gene_bp,
                    "Junction_Reads": row["num_split"],
                    "Spanning_Reads": row["num_paired"]
                }
            )
    return fusions
