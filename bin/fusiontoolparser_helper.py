#!/usr/bin/env python3

"""
Get predicted gene fusions from tools
Fusions are parsed into "id / fusion genes / breakpoints / supporting reads / prediction tool"
and output based on their set recurrency value. The class requires the path' of data and
output folders and writes the "Detected_Fusions.csv"

@author: Tron (PASO)
@version: 20221109
"""
import re
import operator


# TODO: this fixes the reference genome! -> make it more flexible
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
    # Column names as of STAR-fusion v1.12.0
    #FusionName     JunctionReadCount       SpanningFragCount       est_J   est_S   SpliceType      LeftGene        LeftBreakpoint  RightGene       RightBreakpoint JunctionReads   SpanningFrags   LargeAnchorSupport      FFPM    LeftBreakDinuc  LeftBreakEntropy
    #    RightBreakDinuc RightBreakEntropy       annots
    # ETV6--NTRK3     160     390     160.00  390.00  ONLY_REF_SPLICE ETV6^ENSG00000139083.11 chr12:11869969:+        NTRK3^ENSG00000140538.16        chr15:87940753:-
    fusion_map = {}
    with open(infile, "r") as prediction:
        next(prediction)  # skip header line
        for line in prediction:
            elements = line.rstrip().split("\t")

            fusion_gene = elements[0].replace("--", "_").upper()
            up_gene_bp = elements[7].strip("chr")
            dn_gene_bp = elements[9].strip("chr")
            # check whether fusion gene is not on primary chr
            if (
                up_gene_bp.split(":")[0] not in chr_list
                or dn_gene_bp.split(":")[0] not in chr_list
            ):
                continue

            bpid = up_gene_bp + "_" + dn_gene_bp

            fusion_map[bpid] = [
                fusion_gene,  # fusion_gene
                up_gene_bp,  # up_gene_bp
                dn_gene_bp,  # dn_gene_bp
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


def parse_infusion_results(infile):
    """Load and parse results from infusion"""
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


def parse_arriba_results(infile):
    """Load and parse results from arriba"""

    fusion_map = {}
    with open(infile) as prediction:
        next(prediction)  # skip header line
        for line in prediction:
            elements = line.rstrip().split("\t")
            #gene1  gene2   strand1(gene/fusion)    strand2(gene/fusion)    breakpoint1     breakpoint2     site1   site2   type    split_reads1    split_reads2    discordant_mates        coverage1       coverage2       confidence      reading_frame   tags    retained_protein_domains        closest_genomic_breakpoint1     closest_genomic_breakpoint2     gene_id1        gene_id2        transcript_id1  transcript_id2  direction1      direction2      filters fusion_transcript       peptide_sequence        read_identifiers
            # ETV6    NTRK3   +/+     -/-     12:11869969     15:87940753     CDS/splice-site CDS/splice-site translocation   90      40      300     304     301     high    in-frame        Mitelman        Sterile_alpha_motif_(SAM)/Pointed_domain(100%)|Protein_kinase_domain(100%)      .       .       ENSG00000139083 ENSG00000140538 ENST00000396373 ENST00000394480 downstream      downstream      duplicates(34),low_entropy(5),mismatches(10)
            fusion_gene = (elements[0] + "_" + elements[1]).upper()

            # element[4/5] = chr:bp, [2/3] = strands
            # fusion strand will be preferred, but fallback to gene strand should be possible
            up_gene_strand_split = elements[2].split("/")
            dn_gene_strand_split = elements[3].split("/")
            up_gene_strand = ""
            dn_gene_strand = ""
            if up_gene_strand_split[1] != ".":
                up_gene_strand = up_gene_strand_split[1]
            else:
                up_gene_strand = up_gene_strand_split[0]
            
            
            if dn_gene_strand_split[1] != ".":
                dn_gene_strand = dn_gene_strand_split[1]
            else:
                dn_gene_strand = dn_gene_strand_split[0]
            
            up_gene_id = elements[4] + ":" + up_gene_strand
            dn_gene_id = elements[5] + ":" + dn_gene_strand

            # check whether fusion gene is not on primary chr
            if elements[4].split(":")[0] not in chr_list or elements[5].split(":")[0] not in chr_list:
                continue

            bpid = up_gene_id + "_" + dn_gene_id

            junc_reads = str(int(elements[9]) + int(elements[10]))

            fusion_map[bpid] = [
                fusion_gene,  # fusion_gene
                up_gene_id,  # up_gene_bp
                dn_gene_id,  # dn_gene_bp
                junc_reads,  # junc_reads_num
                elements[11],  # span_reads_num
            ]
    return fusion_map
