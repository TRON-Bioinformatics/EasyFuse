"""
Get predicted gene fusions from tools
Fusions are parsed into "id / fusion genes / breakpoints / supporting reads / prediction tool"
and output based on their set recurrency value. The class requires the path' of data and
output folders and writes the "Detected_Fusions.csv"

@author: Tron (PASO)
@version: 20240208
"""

import csv
import operator
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

def parse_fusioncatcher_results(infile_1: str, infile_2: str) -> dict:
    """
    Load and parse results from fusioncatcher
    https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md#6---usage
    """

    reciprocal_fusions = get_reciprocal_fusions(infile_1)
    fusion_map = {}
    # final-list_candidate-fusion-genes.txt
    with open(infile_2, encoding="utf8") as prediction:
        reader = csv.reader(prediction, delimiter="\t")
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

            fusion_map[bpid] = [
                fusion_gene,
                bp1,
                bp2,
                row["Spanning_unique_reads"],
                row["Spanning_pairs"]
            ]
    return fusion_map


def parse_starfusion_results(infile):
    """Load and parse results from star-fusion."""
    fusion_map = {}
    with open(infile, "r", encoding="utf8") as prediction:
        reader = csv.DictReader(prediction, delimiter="\t")
        for row in reader:
            fusion_gene = row["FusionName"].replace("--", "_").upper()
            up_gene_bp = row["LeftBreakpoint"].strip("chr")
            dn_gene_bp = row["RightBreakpoint"].strip("chr")
            # check whether fusion gene is not on primary chr
            if (
                up_gene_bp.split(":")[0] not in CHROMOSOMES
                or dn_gene_bp.split(":")[0] not in CHROMOSOMES
            ):
                continue

            bpid =  + f"{up_gene_bp}_{dn_gene_bp}"

            fusion_map[bpid] = [
                fusion_gene,
                up_gene_bp,
                dn_gene_bp,
                row["JunctionReadCount"],
                row["SpanningFragCount"]
            ]
    return fusion_map


def parse_mapsplice_results(infile):
    """Load and parse results from mapsplice"""
    fusion_map = {}
    with open(infile, encoding="utf8") as prediction:
        for line in prediction:
            elements = line.rstrip().split("\t")

            fusion_gene = (
                elements[60].split(",")[0] + "_" + elements[61].split(",")[0]
            ).upper()

            # element[0] = chr num, [1/2] = breakpoints, [5] = strands
            up_gene_id = (
                elements[0].split("~")[0] + ":" + elements[1] + ":" + elements[5][0]
            )
            dn_gene_id = (
                elements[0].split("~")[1] + ":" + elements[2] + ":" + elements[5][1]
            )

            if (
                up_gene_id.split(":")[0] not in CHROMOSOMES
                or dn_gene_id.split(":")[0] not in CHROMOSOMES
            ):
                continue

            bpid = up_gene_id + "_" + dn_gene_id

            fusion_map[bpid] = [
                fusion_gene,  # fusion_gene
                up_gene_id,  # up_gene_bp
                dn_gene_id,  # dn_gene_bp
                elements[4],  # junc_reads_num
                elements[27],  # span_reads_num
            ]
    return fusion_map


def parse_infusion_results(infile: str) -> dict:
    """
    Load and parse results from infusion
    Description of output here: https://bitbucket.org/kokonech/infusion/wiki/Home
    """
    fusion_map = {}
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

            fusion_map[bpid] = [
                fusion_gene,  # fusion_gene
                up_gene_bp,  # up_gene_bp
                dn_gene_bp,  # dn_gene_bp
                row["num_split"],  # junc_reads_num
                row["num_paired"],  # span_reads_num
            ]
    return fusion_map


def get_fusion_gene_id_infusion(gene_id, transcript_list_field):
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


def parse_soapfuse_results(infile):
    """Load and parse results from soapfuse"""

    fusion_map = {}
    with open(infile, encoding="utf8") as prediction:
        next(prediction)  # skip header line
        for line in prediction:
            elements = line.rstrip().split("\t")

            fusion_gene = (elements[0] + "_" + elements[5]).upper()

            # element[1/6] = chr num, [3/8] = breakpoints, [2/7] = strands
            up_gene_id = elements[1] + ":" + elements[3] + ":" + elements[2]
            dn_gene_id = elements[6] + ":" + elements[8] + ":" + elements[7]

            # check whether fusion gene is not on primary chr
            if elements[1] not in CHROMOSOMES or elements[6] not in CHROMOSOMES:
                continue

            bpid = up_gene_id + "_" + dn_gene_id

            fusion_map[bpid] = [
                fusion_gene,  # fusion_gene
                up_gene_id,  # up_gene_bp
                dn_gene_id,  # dn_gene_bp
                elements[11],  # junc_reads_num
                elements[10],  # span_reads_num
            ]
    return fusion_map


def get_strand(strand_field: str) -> str:
    """
    Helper method for parsing arriba data.
    # strands = gene strand / fusion strand
    # fusion strand will be preferred, but fallback to gene strand should be possible
    """
    strand_split = strand_field.split("/")
    if strand_split[1] != ".":
        return strand_split[1]
    return strand_split[0]


def parse_arriba_results(infile):
    """
    Load and parse results from arriba.
    https://arriba.readthedocs.io/en/latest/output-files/#fusionstsv
    """

    fusion_map = {}
    with open(infile, encoding="utf8") as prediction:
        reader = csv.DictReader(prediction, delimiter="\t")
        for row in reader:
            if row["confidence"] != "high":
                continue
            gene_1 = row["gene1"]
            gene_2 = row["gene2"]
            fusion_gene = f"{gene_1}_{gene_2}".upper()

            up_gene_strand = get_strand(row["strand1(gene/fusion)"])
            dn_gene_strand = get_strand(row["strand2(gene/fusion)"])
            up_gene_bp = row["breakpoint1"] + ":" + up_gene_strand
            dn_gene_bp = row["breakpoint2"] + ":" + dn_gene_strand
            chrom1 = row["breakpoint1"].split(":")[0]
            chrom2 = row["breakpoint2"].split(":")[0]
            # check whether fusion gene is not on primary chr
            if chrom1 not in CHROMOSOMES or chrom2 not in CHROMOSOMES:
                continue

            bpid = f"{up_gene_bp}_{dn_gene_bp}"

            junc_reads = str(int(row["split_reads1"]) + int(row["split_reads2"]))

            fusion_map[bpid] = [
                fusion_gene,
                up_gene_bp,
                dn_gene_bp,
                junc_reads,
                row["discordant_mates"]
            ]
    return fusion_map
