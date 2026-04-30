"""
Module for parsing starfusion results
"""
import csv

from .file_headers import CHROMOSOMES

def parse_starfusion_results(infile: str) -> list:
    """Load and parse results from star-fusion."""
    fusions = []
    with open(infile, "r", encoding="utf8") as prediction:
        reader = csv.DictReader(prediction, delimiter="\t")
        for row in reader:
            fusion_gene = row["#FusionName"].replace("--", "_").upper()
            up_gene_bp = row["LeftBreakpoint"].strip("chr")
            dn_gene_bp = row["RightBreakpoint"].strip("chr")
            # check whether fusion gene is not on primary chr
            if (
                up_gene_bp.split(":")[0] not in CHROMOSOMES
                or dn_gene_bp.split(":")[0] not in CHROMOSOMES
            ):
                continue

            bpid = f"{up_gene_bp}_{dn_gene_bp}"

            fusions.append(
                {
                    "BPID": bpid,
                    "Fusion_Gene": fusion_gene,
                    "Breakpoint1": up_gene_bp,
                    "Breakpoint2": dn_gene_bp,
                    "Junction_Reads": row["JunctionReadCount"],
                    "Spanning_Reads": row["SpanningFragCount"]
                }
            )
    return fusions
