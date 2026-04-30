"""
Module for parsing mapsplice results
"""

import csv

from .file_headers import CHROMOSOMES

def parse_mapsplice_results(infile: str) -> list:
    """Load and parse results from mapsplice"""
    fusions = []
    with open(infile, encoding="utf8") as prediction:
        reader = csv.reader(prediction, delimiter="\t")
        for row in reader:

            fusion_gene = (
                row[60].split(",")[0] + "_" + row[61].split(",")[0]
            ).upper()

            # element[0] = chr num, [1/2] = breakpoints, [5] = strands
            up_gene_id = (
                row[0].split("~")[0] + ":" + row[1] + ":" + row[5][0]
            )
            dn_gene_id = (
                row[0].split("~")[1] + ":" + row[2] + ":" + row[5][1]
            )

            if (
                up_gene_id.split(":")[0] not in CHROMOSOMES
                or dn_gene_id.split(":")[0] not in CHROMOSOMES
            ):
                continue

            bpid = up_gene_id + "_" + dn_gene_id

            fusions.append(
                {
                    "BPID": bpid,
                    "Fusion_Gene": fusion_gene,
                    "Breakpoint1": up_gene_id,
                    "Breakpoint2": dn_gene_id,
                    "Junction_Reads": row[4],
                    "Spanning_Reads": row[27],
                }
            )
    return fusions
