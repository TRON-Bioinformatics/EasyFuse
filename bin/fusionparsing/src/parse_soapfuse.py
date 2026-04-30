"""
Module for parsing soapfuse results
"""

import csv

from .file_headers import CHROMOSOMES

def parse_soapfuse_results(infile: str) -> list:
    """Load and parse results from soapfuse"""

    fusions = []
    with open(infile, encoding="utf8") as prediction:
        reader = csv.DictReader(prediction, delimiter="\t")
        for row in reader:

            fusion_gene = (row["up_gene"] + "_" + row["dw_gene"]).upper()

            # element[1/6] = chr num, [3/8] = breakpoints, [2/7] = strands
            up_gene_id = row["up_chr"] + ":" + row["up_Genome_pos"] + ":" + row["up_strand"]
            dn_gene_id = row["dw_chr"] + ":" + row["dw_Genome_pos"] + ":" + row["dw_strand"]

            # check whether fusion gene is not on primary chr
            if row["up_chr"] not in CHROMOSOMES or row["dw_chr"] not in CHROMOSOMES:
                continue

            bpid = up_gene_id + "_" + dn_gene_id

            fusions.append(
                {
                    "BPID": bpid,
                    "Fusion_Gene": fusion_gene,
                    "Breakpoint1": up_gene_id,
                    "Breakpoint2": dn_gene_id,
                    "Junction_Reads": row["Junc_reads_num"],
                    "Spanning_Reads": row["Span_reads_num"],
                }
            )
    return fusions
