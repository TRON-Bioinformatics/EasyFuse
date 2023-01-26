#!/usr/bin/env python3

"""
Get predicted gene fusions from tools
Fusions are parsed into "id / fusion genes / breakpoints / supporting reads / prediction tool"
and output based on their set recurrency value. The class requires the path' of data and
output folders and writes the "Detected_Fusions.csv"

@author: Tron (PASO)
@version: 20221109
"""
import os
import re
import operator


chr_list = (
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "X",
    "Y",
    "MT",
)


def parse_results(output_path, tool):
    """Load and parse results from tool"""

    fusioncatcher_results_path_1 = os.path.join(
        output_path, "fusion", "fusioncatcher", "summary_candidate_fusions.txt"
    )
    fusioncatcher_results_path_2 = os.path.join(
        output_path, "fusion", "fusioncatcher", "final-list_candidate-fusion-genes.txt"
    )
    starfusion_results_path = os.path.join(
        output_path, "fusion", "starfusion", "star-fusion.fusion_predictions.tsv"
    )
    mapsplice_results_path = os.path.join(
        output_path, "fusion", "mapsplice", "fusions_well_annotated.txt"
    )
    infusion_results_path = os.path.join(
        output_path, "fusion", "infusion", "fusions.detailed.txt"
    )
    soapfuse_results_path = ""
    folder_to_scan = os.path.join(
        output_path, "fusion", "soapfuse", "final_fusion_genes"
    )
    for filename in os.listdir(folder_to_scan):
        folder_path = os.path.join(folder_to_scan, filename)
        if os.path.isdir(folder_path):
            for res in os.listdir(folder_path):
                if res.endswith(".final.Fusion.specific.for.genes"):
                    soapfuse_results_path = os.path.join(folder_path, res)

    starchip_results_path = os.path.join(
        output_path, "fusion", "starchip", "starchip.summary"
    )

    if tool == "fusioncatcher":
        return parse_fusioncatcher_results(
            fusioncatcher_results_path_1, fusioncatcher_results_path_2
        )
    elif tool == "starfusion":
        return parse_starfusion_results(starfusion_results_path)
    elif tool == "mapsplice":
        return parse_mapsplice_results(mapsplice_results_path)
    elif tool == "starchip":
        return parse_starchip_results(starchip_results_path)
    elif tool == "infusion":
        return parse_infusion_results(infusion_results_path)
    elif tool == "soapfuse":
        return parse_soapfuse_results(soapfuse_results_path)
    return None


def parse_fusioncatcher_results(infile1, infile2):
    """Load and parse results from fusioncatcher"""

    reciprocal_fusions = []
    with open(infile1) as predict_summary:
        for line in predict_summary:
            if line.strip().startswith("*"):
                fusion_gene = re.search(r"\*\s([\S]*)", line).group(1)
                if "reciprocal" in line:
                    reciprocal_fusions.append(fusion_gene.replace("--", "_").upper())
    fusion_map = {}
    with open(infile2) as prediction:
        next(prediction)  # skip header line
        for line in prediction:
            elements = line.rstrip().split("\t")

            fusion_gene = (elements[0] + "_" + elements[1]).upper()
            # if the fusion gene is reciprocal, the fusion id is reversed? <- what for??
            if fusion_gene in reciprocal_fusions:
                fusion_gene = (elements[1] + "_" + elements[0]).upper()

            # skip all prediction not on standard chromosomes
            if (
                elements[8].split(":")[0] not in chr_list
                or elements[9].split(":")[0] not in chr_list
            ):
                continue

            bpid = elements[8] + "_" + elements[9]

            fusion_map[bpid] = [
                fusion_gene,  # fusion_gene
                elements[8],  # up_gene_bp
                elements[9],  # dn_gene_bp
                elements[5],  # junc_reads_num
                elements[4],  # span_reads_num
            ]
    return fusion_map


def parse_starfusion_results(infile):
    """Load and parse results from star-fusion"""
    fusion_map = {}
    with open(infile, "r") as prediction:
        next(prediction)  # skip header line
        for line in prediction:
            elements = line.rstrip().split("\t")

            fusion_gene = elements[0].replace("--", "_").upper()
            # check whether fusion gene is not on primary chr
            if (
                elements[5].split(":")[0] not in chr_list
                or elements[7].split(":")[0] not in chr_list
            ):
                continue

            bpid = elements[5] + "_" + elements[7]

            fusion_map[bpid] = [
                fusion_gene,  # fusion_gene
                elements[5],  # up_gene_bp
                elements[7],  # dn_gene_bp
                elements[1],  # junc_reads_num
                elements[2],  # span_reads_num
            ]
    return fusion_map


def parse_mapsplice_results(infile):
    """Load and parse results from mapsplice"""
    fusion_map = {}
    with open(infile) as prediction:
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
                up_gene_id.split(":")[0] not in chr_list
                or dn_gene_id.split(":")[0] not in chr_list
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


def parse_starchip_results(infile):
    """Load and parse results from starchip"""
    fusion_map = {}
    with open(infile) as prediction:
        next(prediction)  # skip header line
        for line in prediction:
            elements = line.rstrip().split("\t")

            fusion_gene = "{0}_{1}".format(elements[5], elements[7]).upper()

            # check whether fusion gene is not on primary chr
            if (
                elements[0].split(":")[0] not in self.chr_list
                or elements[1].split(":")[0] not in self.chr_list
            ):
                continue

            bpid = elements[0] + "_" + elements[1]

            fusion_map[bpid] = [
                fusion_gene,  # fusion_gene
                elements[0],  # up_gene_bp
                elements[1],  # dn_gene_bp
                elements[3],  # junc_reads_num
                elements[2],  # span_reads_num
                "Starchip",
            ]
    return fusion_map


def parse_infusion_results(infile):
    """Load and parse results from starchip"""
    fusion_map = {}
    with open(infile) as prediction:
        next(prediction)  # skip header line
        for line in prediction:
            elements = line.rstrip().split("\t")

            fusion_gene = "{0}_{1}".format(
                get_fusion_gene_id_infusion(elements[21], elements[22]),
                get_fusion_gene_id_infusion(elements[27], elements[28]),
            ).upper()

            # element[1/4] = chr num, [2/5] = breakpoints, [23/29] = strands
            up_gene_id = elements[1] + ":" + elements[2] + ":" + elements[23]
            dn_gene_id = elements[4] + ":" + elements[5] + ":" + elements[29]

            # check whether fusion gene is not on primary chr
            if elements[1] not in chr_list or elements[4] not in chr_list:
                continue

            bpid = up_gene_id + "_" + dn_gene_id

            fusion_map[bpid] = [
                fusion_gene,  # fusion_gene
                up_gene_id,  # up_gene_bp
                dn_gene_id,  # dn_gene_bp
                elements[7],  # junc_reads_num
                elements[8],  # span_reads_num
            ]
    return fusion_map


def get_fusion_gene_id_infusion(gene_id, transcript_list_field):
    """Helper method for infusion data parsing. Returns the most frequently listed fusion id in a list of fusion ids"""
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
    with open(infile) as prediction:
        next(prediction)  # skip header line
        for line in prediction:
            elements = line.rstrip().split("\t")

            fusion_gene = (elements[0] + "_" + elements[5]).upper()

            # element[1/6] = chr num, [3/8] = breakpoints, [2/7] = strands
            up_gene_id = elements[1] + ":" + elements[3] + ":" + elements[2]
            dn_gene_id = elements[6] + ":" + elements[8] + ":" + elements[7]

            # check whether fusion gene is not on primary chr
            if elements[1] not in chr_list or elements[6] not in chr_list:
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
