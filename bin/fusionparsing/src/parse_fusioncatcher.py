"""
Module for parsing fusioncatcher results
"""

import csv
import re

from .file_headers import CHROMOSOMES


def get_reciprocal_fusions(infile: str) -> list:
    """
    Load and parse results from fusioncatcher:
    summary_candidate_fusions.txt
    """
    reciprocal_fusions = []
    with open(infile, encoding="utf8") as predict_summary:
        for line in predict_summary:
            if line.strip().startswith("*"):
                fusion_gene = re.search(r"\*\s([\S]*)", line).group(1)
                if "reciprocal" in line:
                    reciprocal_fusions.append(fusion_gene.replace("--", "_").upper())
    return reciprocal_fusions


def parse_fusioncatcher_results(infile_1: str, infile_2: str) -> list:
    """
    Load and parse results from fusioncatcher
    https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md#6---usage
    """

    reciprocal_fusions = get_reciprocal_fusions(infile_1)
    fusions = []
    # final-list_candidate-fusion-genes.txt
    with open(infile_2, encoding="utf8") as prediction:
        reader = csv.DictReader(prediction, delimiter="\t")
        for row in reader:
            gene_1 = row["Gene_1_symbol(5end_fusion_partner)"]
            gene_2 = row["Gene_2_symbol(3end_fusion_partner)"]
            fusion_gene = f"{gene_1}_{gene_2}".upper()
            # if the fusion gene is reciprocal, the fusion id is reversed? <- what for??
            if fusion_gene in reciprocal_fusions:
                fusion_gene = f"{gene_2}_{gene_1}".upper()

            bp1 = row["Fusion_point_for_gene_1(5end_fusion_partner)"]
            bp2 = row["Fusion_point_for_gene_2(3end_fusion_partner)"]
            # skip all prediction not on standard chromosomes
            if (
                bp1.split(":")[0] not in CHROMOSOMES
                or bp2.split(":")[0] not in CHROMOSOMES
            ):
                continue

            bpid = f"{bp1}_{bp2}"

            fusions.append(
                {
                    "BPID": bpid,
                    "Fusion_Gene": fusion_gene,
                    "Breakpoint1": bp1,
                    "Breakpoint2": bp2,
                    "Junction_Reads": row["Spanning_unique_reads"],
                    "Spanning_Reads": row["Spanning_pairs"],
                }
            )
    return fusions
