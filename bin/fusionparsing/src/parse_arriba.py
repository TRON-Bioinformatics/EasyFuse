"""
Module for parsing arriba results
"""

import csv

from .file_headers import CHROMOSOMES


def parse_strand(strand_field: str) -> str:
    """
    Helper method for parsing arriba data.
    # strands = gene strand / fusion strand
    # fusion strand will be preferred, but fallback to gene strand should be possible
    """
    strand_split = strand_field.split("/")
    if strand_split[1] != ".":
        return strand_split[1]
    return strand_split[0]


def parse_arriba_results(infile: str) -> list:
    """
    Load and parse results from arriba.
    https://arriba.readthedocs.io/en/latest/output-files/#fusionstsv
    """

    fusions = []
    with open(infile, encoding="utf8") as prediction:
        reader = csv.DictReader(prediction, delimiter="\t")
        for row in reader:
            if row["confidence"] != "high":
                continue
            gene_1 = row["#gene1"]
            gene_2 = row["gene2"]
            fusion_gene = f"{gene_1}_{gene_2}".upper()

            up_gene_strand = parse_strand(row["strand1(gene/fusion)"])
            dn_gene_strand = parse_strand(row["strand2(gene/fusion)"])
            up_gene_bp = row["breakpoint1"] + ":" + up_gene_strand
            dn_gene_bp = row["breakpoint2"] + ":" + dn_gene_strand
            chrom1 = row["breakpoint1"].split(":")[0]
            chrom2 = row["breakpoint2"].split(":")[0]
            # check whether fusion gene is not on primary chr
            if chrom1 not in CHROMOSOMES or chrom2 not in CHROMOSOMES:
                continue

            bpid = f"{up_gene_bp}_{dn_gene_bp}"

            junc_reads = str(int(row["split_reads1"]) + int(row["split_reads2"]))

            fusions.append(
                {
                    "BPID": bpid,
                    "Fusion_Gene": fusion_gene,
                    "Breakpoint1": up_gene_bp,
                    "Breakpoint2": dn_gene_bp,
                    "Junction_Reads": junc_reads,
                    "Spanning_Reads": row["discordant_mates"]
                }
            )
    return fusions
