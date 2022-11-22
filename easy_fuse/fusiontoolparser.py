#!/usr/bin/env python3

"""
Get predicted gene fusions from tools
Fusions are parsed into "id / fusion genes / breakpoints / supporting reads / prediction tool"
and output based on their set recurrency value. The class requires the path' of data and
output folders and writes the "Detected_Fusions.csv"

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""

import operator
import os
import re
import sys
from argparse import ArgumentParser

import logzero
from logzero import logger


# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class FusionParser(object):
    """Get and parse results from previously run programs (fusion prediction, hla typing, expression estimation)"""

    # Initialization of parameters
    # urla: todo: pylint convention - too many arguments
    #             are all vars required to be class vars?
    def __init__(self, scratch_path, fusion_output_path, sample_id, tool_num_cutoff, fusiontool_list, sample_log):
        """Parameter initiation and work folder creation."""
        self.scratch_path = scratch_path
        self.fusion_output_path = fusion_output_path
        self.sample_id = sample_id
        self.tool_num_cutoff = int(tool_num_cutoff)
        # urla: if we want to be more generic and allow different annotations, identification of the chr names
        #       (eg "chr1" vs "1" and "chrM" vs "MT") should be performed in advance
        self.chr_list = (
            "1", "2", "3", "4", "5",
            "6", "7", "8", "9", "10",
            "11", "12", "13", "14", "15",
            "16", "17", "18", "19", "20",
            "21", "22", "X", "Y", "MT"
        )
        self.tools = fusiontool_list.split(",")
        logzero.logfile(sample_log)

    #####
    ### fusion tool parser
    #
    # fusioncatcher - results file is "summary_candidate_fusions.txt"
    def get_fusioncatcher_results(self):
        """Load and parse results from fusioncatcher"""
        fusioncatcher_predict_summary = os.path.join(self.scratch_path, "fusion", "fusioncatcher",
                                                     "summary_candidate_fusions.txt")
        fusioncatcher_predict_detail = os.path.join(self.scratch_path, "fusion", "fusioncatcher",
                                                    "final-list_candidate-fusion-genes.txt")
        reciprocal_fusions = []
        with open(fusioncatcher_predict_summary) as predict_summary:
            for line in predict_summary:
                if line.strip().startswith("*"):
                    # urla: todo: pylint warning - anomalous backslash
                    #             I don't really understand the problem and the RE is simple and working very fine
                    match = re.search(r'\*\s([\S]*)', line)
                    if match:
                        fusion_gene = match.group(1)
                        if "reciprocal" in line:
                            reciprocal_fusions.append(fusion_gene.replace("--", "_").upper())
        fusion_map = {}
        with open(fusioncatcher_predict_detail) as prediction:
            next(prediction)  # skip header line
            for line in prediction:
                elements = line.rstrip().split("\t")
                # Currently relevant fields (marked *) in the output file are:
                # * elements[0] ~ Gene_1_symbol(5end_fusion_partner)
                # * elements[1] ~ Gene_2_symbol(3end_fusion_partner)
                #   elements[2] ~ Fusion_description
                #   elements[3] ~ Counts_of_common_mapping_reads
                # * elements[4] ~ Spanning_pairs
                # * elements[5] ~ Spanning_unique_reads
                #   elements[6] ~ Longest_anchor_found
                #   elements[7] ~ Fusion_finding_method
                # * elements[8] ~ Fusion_point_for_gene_1(5end_fusion_partner)
                # * elements[9] ~ Fusion_point_for_gene_2(3end_fusion_partner)
                #   elements[10] ~ Gene_1_id(5end_fusion_partner)
                #   elements[11] ~ Gene_2_id(3end_fusion_partner)
                #   elements[12] ~ Exon_1_id(5end_fusion_partner)
                #   elements[13] ~ Exon_2_id(3end_fusion_partner)
                #   elements[14] ~ Fusion_sequence
                #   elements[15] ~ Predicted_effect
                fusion_gene = (elements[0] + "_" + elements[1]).upper()
                # if the fusion gene is reciprocal, the fusion id is reversed? <- what for??
                if fusion_gene in reciprocal_fusions:
                    fusion_gene = (elements[1] + "_" + elements[0]).upper()
                #                for key in self.sub_dict:
                #                    fusion_gene = fusion_gene.replace(key, self.sub_dict[key])

                # urla: why not catching "Counts_of_common_mapping_reads" which indicate how similar the fusion partners are? <- FP's; should be 0 for max specificity
                # common_map_num = elements[3]
                # urla: according to online manual, junction reads are in elements[5] and all supporting (junction and spanning?) in elements[4]?!
                # urla_c: I changed it here

                # skip all prediction not on standard chromosomes
                if elements[8].split(":")[0] not in self.chr_list or elements[9].split(":")[0] not in self.chr_list:
                    continue

                bpid = elements[8] + "_" + elements[9]

                fusion_map[bpid] = [
                    fusion_gene,  # fusion_gene
                    elements[8],  # up_gene_bp
                    elements[9],  # dn_gene_bp
                    elements[5],  # junc_reads_num - urla: todo: verify that this is correct
                    elements[4],  # span_reads_num - urla: todo: verify that this is correct
                    self.sample_id,
                    "Fusioncatcher"
                ]
        return fusion_map

    # starfusion - results file is "star-fusion.fusion_candidates.final (not any more, now: star-fusion.fusion_predictions.abridged.tsv)"
    # in order to be compatible with old and new starfusion versions, the new file name is checked first, if not available, the old one is used.
    # urla - note: this is rather a hack then a proper solution, but proper version handling would actually require to dig through all versions
    #              and look for changes in file names and output columns and maybe additional stuff... this would properly require a tremendous
    #              amount of time and is therefore out of scope of this project
    def get_starfusion_results(self):
        """Load and parse results from star-fusion"""
        starfusion_predict = os.path.join(self.scratch_path, "fusion", "starfusion",
                                          "star-fusion.fusion_predictions.abridged.tsv")
        if not os.path.isfile(starfusion_predict):
            starfusion_predict = os.path.join(self.scratch_path, "fusion", "starfusion",
                                              "star-fusion.fusion_candidates.final")
        fusion_map = {}
        with open(starfusion_predict, "r") as prediction:
            next(prediction)  # skip header line
            for line in prediction:
                elements = line.rstrip().split("\t")
                # Currently relevant fields (marked *) in the output file are:
                # * elements[0] ~ FusionName
                # * elements[1] ~ JunctionReadCount
                # * elements[2] ~ SpanningFragCount
                #   elements[3] ~ SpliceType
                #   elements[4] ~ LeftGene
                # * elements[5] ~ LeftBreakpoint
                #   elements[6] ~ RightGene
                # * elements[7] ~ RightBreakpoint
                #   elements[8] ~ LargeAnchorSupport
                #   elements[9] ~ FFPM
                #   elements[10] ~ LeftBreakDinuc
                #   elements[11] ~ LeftBreakEntropy
                #   elements[12] ~ RightBreakDinuc
                #   elements[13] ~ RightBreakEntropy
                #   elements[14] ~ annots
                fusion_gene = elements[0].replace("--", "_").upper()
                # check whether fusion gene is not on primary chr
                if elements[5].split(":")[0] not in self.chr_list or elements[7].split(":")[0] not in self.chr_list:
                    continue
                bpid = elements[5] + "_" + elements[7]

                fusion_map[bpid] = [
                    fusion_gene,  # fusion_gene
                    elements[5],  # up_gene_bp
                    elements[7],  # dn_gene_bp
                    elements[1],  # junc_reads_num
                    elements[2],  # span_reads_num
                    self.sample_id,
                    "Starfusion"
                ]
        return fusion_map

    # mapsplice2 - results file is "fusions_well_annotated.txt"
    def get_mapsplice_results(self):
        """Load and parse results from mapsplice"""
        mapsplice_predict = os.path.join(self.scratch_path, "fusion", "mapsplice", "fusions_well_annotated.txt")
        fusion_map = {}
        with open(mapsplice_predict) as prediction:
            # next(prediction) # mapsplice final result table, doesn't have a header!
            for line in prediction:
                elements = line.rstrip().split("\t")
                # Currently relevant fields (marked *) in the output file are:
                # * elements[0] ~ chrom: the two chromosomes involved in fusion junction
                # * elements[1] ~ doner_end: The end position of doner site of splicing on chromosome
                # * elements[2] ~ acceptor_start: The start position of acceptor site of splicing on chromosome
                #   elements[3] ~ id: The id of fusion junction
                # * elements[4] ~ coverage: number of reads aligned to the fusion junction
                # * elements[5] ~ strand: strand of the reads mapped to the two chromosomes
                #   elements[6] ~ rgb: An RGB value of the form R,G,B
                #   elements[7] ~ block_count:  The number of blocks in the BED line
                #   elements[8] ~ block_size: A comma-separated list of the block sizes.
                #   elements[9] ~ block_distance: A comma-separated list of block distance.
                #   elements[10] ~ entropy: entropy of the fusion junction.
                #   elements[11] ~ flank_case: non-zero for canonical and semi-canonical junctions (ATAC 1;GTAT 2;CTGC 3;GCAG 4;GTAG 5;CTAC 6;others 0)
                #   elements[12] ~ flank_string: the two basepairs after doner site combined the two basepairs before acceptor site
                #   elements[13] ~ min_mismatch: Minimal mismatch of read mapped to the fusion junction
                #   elements[14] ~ max_mismatch: Maximal mismatch of read mapped to the fusion junction
                #   elements[15] ~ ave_mismatch: Average mismatch of all reads mapped to the junction
                #   elements[16] ~ max_min_suffix:  if doner site is shorter than acceptor site, and if the doner site is longer than current maximal doner site length, then update current maximal doner site length
                #   elements[17] ~ max_min_prefix: if acceptor site is shorter than doner site, and if the doner site is longer than current maximal acceptor site length, then update current maximal acceptor site length
                #   elements[18] ~ min_anchor_difference: Minimal difference between doner site and acceptor site
                # * elements[19] ~ unique_read_count: Number of uniquely mapped reads mapped to the fusion
                #   elements[20] ~ multi_read_count: Number of multiple mapped reads mapped to the fusion
                #   elements[21] ~ paired_read_count: Number of reads mapped to fusion and can be paired with their mates near the fusion
                #   elements[22] ~ left_paired_read_count: Number of paired reads that the read itself is mapped to the left of its mate on genome
                #   elements[23] ~ right_paired_read_count: Number of paired reads that the read itself is mapped to the right of its mate on genome
                #   elements[24] ~ multiple_paired_read_count: Number of multiple mapped reads mapped to the fusion and are paired with their mates
                #   elements[25] ~ unique_paired_read_count: Number of uniquely mapped reads mapped to the fusion and are paired with their mates
                #   elements[26] ~ single_read_count: Number of reads mapped to the fusion but can't be paired with their mates
                #   elements[27] ~ encompassing_read pair_count: Number of reads pairs surround the fusion(but not cross the fusion)
                #   elements[28] ~ doner_start: The start of doner site of splicing on chromosome
                #   elements[29] ~ acceptor_end: The end of acceptor site of splicing on chromosome
                #   elements[30] ~ doner_iosforms: The isoform(transcript) structure on the doner site. each isoform structure is separated by '|'. The format of each isoform structure is the "start_of_the_isoform,CIGAR_string_of_structure". E.g. 59445681,180M12006N66M8046N47M|59445681,180M20118N47M| Two isoforms start at 59445681
                #   elements[31] ~ acceptor_isoforms: The isoform(transcript) structure on the acceptor site.
                #   elements[32] ~ doner uniformity score (obsolete): The p-value of T-test against the hypothesis that the start of spanning read pairs and encompassing read pairs distribute uniformly on the doner site
                #   elements[33] ~ acceptor uniformity score (obsolete): The p-value of Kolmogorov-Smirnov test against the hypothesis that the end of spanning read pairs and encompassing read pairs distribute uniformly on the acceptor site
                #   elements[34] ~ doner uniformity KS-test score (obsolete): The score of Kolmogorov-Smirnov test against the hypothesis that the start of spanning read pairs and encompassing read pairs distribute uniformly on the doner site
                #   elements[35] ~ acceptor uniformity KS-test score (obsolete): The score of Kolmogorov-Smirnov test against the hypothesis that the end of spanning read pairs and encompassing read pairs distribute uniformly on the acceptor site
                #   elements[36] ~ minimal_doner_isoform_length: Minimal length of isoform structure on the doner site
                #   elements[37] ~ maximal_doner_isoform_length: Maximal length of isoform structure on the doner site
                #   elements[38] ~ minimal_acceptor_isoform_length: Minimal length of isoform structure on the acceptor site
                #   elements[39] ~ maximal_acceptor_isoform_length: Maximal length of isoform structure on the acceptor site
                #   elements[40] ~ paired_reads_entropy: entropy of different read pairs
                #   elements[41] ~ mismatch_per_bp: Average mismatch per base.
                #   elements[42] ~ anchor_score: Anchor score.
                #   elements[43] ~ max_doner_fragment: Maximum doner fragment size.
                #   elements[44] ~ max_acceptor_fragment: Maximum acceptor fragment size.
                #   elements[45] ~ max_cur_fragment: Maximum total fragment length of doner and acceptor.
                #   elements[46] ~ min_cur_fragment: Minimum total fragment length of doner and acceptor
                #   elements[47] ~ ave_cur_fragment: Average total fragment length of doner and acceptor.
                #   elements[48] ~ doner_encompass_unique: Number of uniquely mapped reads surround the donor site of the fusion
                #   elements[49] ~ doner_encompass_multiple: Number of multiple mapped reads surround the donor site of the fusion
                #   elements[50] ~ acceptor_encompass_unique: Number of uniquely mapped reads surround the acceptor site of the fusion
                #   elements[51] ~ acceptor_encompass_multiple: Number of multiple mapped reads surround the acceptor site of the fusion
                #   elements[52] ~ doner_match_to_normal: If the fusion doner site is matched to a normal splice junction. 1 is matched, 0 is not matched
                #   elements[53] ~ acceptor_match_to_normal: If the fusion doner site is matched to a normal splice junction. 1 is matched, 0 is not matched
                #   elements[54] ~ doner_seq: The 25bp sequence at doner site matched to the fusion reads, doner end base included.
                #                  if doner strand is +, it is chrom1[donerEnd-24:donerEnd]
                #                  if doner strand  is -, it is revcomp(chrom1[donerEnd:donerEnd+24]))
                #   elements[55] ~ acceptor_seq: The 25bp sequence at acceptor site matched to the fusion reads, acceptor start base included.
                #                  if acceptor strand is +, it is chrom2[acceptorStart:acceptorStart+24]
                #                  if acceptor strand is -, it is revcomp(chrom2[acceptorStart-24:acceptorStart])
                #                  Note: due to our internal use purposes, the acceptor_seq is always reverse complemented again in MapSplice fusion junction file, which makes it effectively as the following. You can do a reverse completement to acceptor_seq to make the sequence exactly the same as it is mentioned above:
                #                  if acceptor strand is +, it is revcom(chrom2[acceptorStart:acceptorStart+24])
                #                  if acceptor strand is -, it is chrom2[acceptorStart-24:acceptorStart]
                #   elements[56] ~ match_gene_strand (only if --gene-gtf specified): If the fusion strand matched with the annotated gene strand. 1 is matched, 0 is not matched
                #   elements[57] ~ annotated_type (only if --gene-gtf specified): The source of fusion
                #                  from_fusion: The fusion is from fusion alignments
                #                  from_normal: The fusion is from normal alignments, which is normal junction cross two genes(read through fusions)
                #   elements[58] ~ fusion_type (only if --gene-gtf specified): The type of fusion based on the annotated gene
                #                  fusion: The start and end of the fusion is annotated to two distinct genes
                #                  normal: The start and end of the fusion is annotated to same gene(Circular RNAs)
                #                  intergenic: Either the start or end has no gene annotated
                #                  overlapped?
                #   elements[59] ~ gene_strand (only if --gene-gtf specified): The annotated genes strands, if there are
                # * elements[60] ~ annotated_gene_donor (only if --gene-gtf specified): The name of the gene annotated to the doner site of the fusion
                # * elements[61] ~ annotated_gene_acceptor (only if --gene-gtf specified): The name of the gene annotated to the acceptor site of the fusion
                fusion_gene = (elements[60].split(",")[0] + "_" + elements[61].split(",")[0]).upper()

                # element[0] = chr num, [1/2] = breakpoints, [5] = strands
                up_gene_id = elements[0].split("~")[0] + ":" + elements[1] + ":" + elements[5][0]
                dn_gene_id = elements[0].split("~")[1] + ":" + elements[2] + ":" + elements[5][1]

                if up_gene_id.split(":")[0] not in self.chr_list or dn_gene_id.split(":")[0] not in self.chr_list:
                    continue

                bpid = up_gene_id + "_" + dn_gene_id

                fusion_map[bpid] = [
                    fusion_gene,  # fusion_gene
                    up_gene_id,  # up_gene_bp
                    dn_gene_id,  # dn_gene_bp
                    elements[4],  # junc_reads_num
                    elements[27],  # span_reads_num
                    self.sample_id,
                    "Mapsplice"
                ]
        return fusion_map

    # starchip - results file is "starchip.summary"
    def get_starchip_results(self):
        """Load and parse results from starchip"""
        starchip_predict = os.path.join(self.scratch_path, "fusion", "starchip", "starchip.summary")
        fusion_map = {}
        with open(starchip_predict, "r") as prediction:
            next(prediction)  # skip header line
            for line in prediction:
                elements = line.rstrip().split("\t")
                # Currently relevant fields (marked *) in the output file are:
                # * elements[0] ~ Partner1
                # * elements[1] ~ Partner2
                # * elements[2] ~ SpanningReads
                # * elements[3] ~ SplitReads
                #   elements[4] ~ AvgAS
                # * elements[5] ~ NearGene1
                #   elements[6] ~ Distance1
                # * elements[7] ~ NearGene2
                #   elements[8] ~ Distance2
                #   elements[9] ~ ConsensusSeq
                fusion_gene = "{0}_{1}".format(elements[5], elements[7]).upper()

                # check whether fusion gene is not on primary chr
                if elements[0].split(":")[0] not in self.chr_list or elements[1].split(":")[0] not in self.chr_list:
                    continue

                bpid = elements[0] + "_" + elements[1]

                fusion_map[bpid] = [
                    fusion_gene,  # fusion_gene
                    elements[0],  # up_gene_bp
                    elements[1],  # dn_gene_bp
                    elements[3],  # junc_reads_num
                    elements[2],  # span_reads_num
                    self.sample_id,
                    "Starchip"
                ]
        return fusion_map

    # infusion - results file is "fusions.detailed.txt"
    def get_infusion_results(self):
        """Load and parse results from starchip"""
        infusion_predict = os.path.join(self.scratch_path, "fusion", "infusion", "fusions.detailed.txt")
        fusion_map = {}
        with open(infusion_predict, "r") as prediction:
            next(prediction)  # skip header line
            for line in prediction:
                elements = line.rstrip().split("\t")
                # Currently relevant fields (marked *) in the output file are:
                #   elements[0] ~ id
                #   elements[1] ~ ref1
                #   elements[2] ~ break_pos1
                #   elements[3] ~ region1
                #   elements[4] ~ ref2
                #   elements[5] ~ break_pos2
                #   elements[6] ~ region2
                #   elements[7] ~ num_split
                #   elements[8] ~ num_paired
                #   elements[9] ~ num_split_with_pair
                #   elements[10] ~ num_split_rescued
                #   elements[11] ~ num_uniq_starts
                #   elements[12] ~ pap_rate
                #   elements[13] ~ mean_split_pos
                #   elements[14] ~ split_pos_std
                #   elements[15] ~ homogeneity
                #   elements[16] ~ coverage_context
                #   elements[17] ~ ssp
                #   elements[18] ~ fusion_class
                #   elements[19] ~ break_on_exon
                #   elements[20] ~ feature_1
                #   elements[21] ~ gene_1
                #   elements[22] ~ transcript_1
                #   elements[23] ~ gene_1_strand
                #   elements[24] ~ biotype_1
                #   elements[25] ~ expression_1
                #   elements[26] ~ feature_2
                #   elements[27] ~ gene_2
                #   elements[28] ~ transcript_2
                #   elements[29] ~ gene_2_strand
                #   elements[30] ~ biotype_2
                #   elements[31] ~ expression_2
                #   elements[32] ~ splice_motif
                #   elements[33] ~ filters
                fusion_gene = "{0}_{1}".format(self.get_fusion_gene_id_infusion(elements[21], elements[22]),
                                               self.get_fusion_gene_id_infusion(elements[27], elements[28])).upper()

                # element[1/4] = chr num, [2/5] = breakpoints, [23/29] = strands
                up_gene_id = elements[1] + ":" + elements[2] + ":" + elements[23]
                dn_gene_id = elements[4] + ":" + elements[5] + ":" + elements[29]

                # check whether fusion gene is not on primary chr
                if elements[1] not in self.chr_list or elements[4] not in self.chr_list:
                    continue

                bpid = up_gene_id + "_" + dn_gene_id

                fusion_map[bpid] = [
                    fusion_gene,  # fusion_gene
                    up_gene_id,  # up_gene_bp
                    dn_gene_id,  # dn_gene_bp
                    elements[7],  # junc_reads_num
                    elements[8],  # span_reads_num
                    self.sample_id,
                    "Infusion"
                ]
        return fusion_map

    @staticmethod
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

        best_hit = max(gid_dict.iteritems(), key=operator.itemgetter(1))[0]

        return best_hit.replace("_", "-")

    # soapfuse - results file is "*.final.Fusion.specific.for.genes" in "final_fusion_genes"
    # urla - note: it seems like soapfuse is prefixing chromsome ids with "chr" even if eg ensembl data is used which does not have this
    #              I'm replacing this here because this is our recommended dataset...
    def get_soapfuse_results(self):
        """Load and parse results from soapfuse"""
        soapfuse_predict = ""
        folder_to_scan = os.path.join(self.scratch_path, "fusion", "soapfuse", "final_fusion_genes")
        for filename in os.listdir(folder_to_scan):
            folder_path = os.path.join(folder_to_scan, filename)
            if os.path.isdir(folder_path):
                for res in os.listdir(folder_path):
                    if res.endswith(".final.Fusion.specific.for.genes"):
                        soapfuse_predict = os.path.join(folder_path, res)
        if not soapfuse_predict:
            soapfuse_predict = os.path.join(self.scratch_path, "fusion", "soapfuse", "final_fusion_genes",
                                            self.sample_id, self.sample_id + ".final.Fusion.specific.for.genes")
        fusion_map = {}
        with open(soapfuse_predict) as prediction:
            next(prediction)  # skip header line
            for line in prediction:
                elements = line.rstrip().split("\t")
                # Currently relevant fields (marked *) in the output file are:
                # * elements[0] ~ up_gene
                # * elements[1] ~ up_chr
                # * elements[2] ~ up_strand
                # * elements[3] ~ up_Genome_pos
                #   elements[4] ~ up_loc
                # * elements[5] ~ dw_gene
                # * elements[6] ~ dw_chr
                # * elements[7] ~ dw_strand
                # * elements[8] ~ dw_Genome_pos
                #   elements[9] ~ dw_loc
                # * elements[10] ~ Span_reads_num
                # * elements[11] ~ Junc_reads_num
                #   elements[12] ~ Fusion_Type
                #   elements[13] ~ down_fusion_part_frame-shift_or_not
                fusion_gene = (elements[0] + "_" + elements[5]).upper()

                # element[1/6] = chr num, [3/8] = breakpoints, [2/7] = strands
                up_gene_id = elements[1] + ":" + elements[3] + ":" + elements[2]
                dn_gene_id = elements[6] + ":" + elements[8] + ":" + elements[7]

                # check whether fusion gene is not on primary chr
                if elements[1] not in self.chr_list or elements[6] not in self.chr_list:
                    continue

                bpid = up_gene_id + "_" + dn_gene_id

                fusion_map[bpid] = [
                    fusion_gene,  # fusion_gene
                    up_gene_id,  # up_gene_bp
                    dn_gene_id,  # dn_gene_bp
                    elements[11],  # junc_reads_num
                    elements[10],  # span_reads_num
                    self.sample_id,
                    "Soapfuse"
                ]
        return fusion_map

    #####
    ### parser of fusion tool parser outputs
    #
    # urla: at least writing to tool_state_path is sort of obsolete as it is never used downstream of here
    def concatenate_fusion_results(self, tool_state_path, fusion_output_path):
        """Return tuple of (dict of results dicts, dict of errors) and writes error/pass to summary file"""
        with open(tool_state_path, "a") as outf:
            fusion_result_dict = {}  # stores the complete fusion map per tool
            results_with_errors_dict = {}  # stores booleans whether exeption was raised during tool exe

            logger.info("Processing " + self.sample_id)
            sample_string = self.sample_id
            for tool in self.tools:
                fusion_result_dict[tool], results_with_errors_dict[tool] = self.get_tool_results(fusion_output_path,
                                                                                                 tool)

                if results_with_errors_dict[tool]:
                    sample_string += ";0"
                else:
                    sample_string += ";1"
            sample_string.rstrip(";")
            outf.write(sample_string + "\n")
        return (fusion_result_dict, results_with_errors_dict)

    def get_tool_results(self, output_folder_path, tool):
        """Return tuple of (dict of results from individual fusion tools, error type)"""
        pred_res_dict = {}  # dictionary of prediction results
        logger.info("Parsing results for " + tool)
        try:
            if tool == "Fusioncatcher":
                pred_res_dict = self.get_fusioncatcher_results()
            elif tool == "Starfusion":
                pred_res_dict = self.get_starfusion_results()
            elif tool == "Mapsplice":
                pred_res_dict = self.get_mapsplice_results()
            elif tool == "Starchip":
                pred_res_dict = self.get_starchip_results()
            elif tool == "Infusion":
                pred_res_dict = self.get_infusion_results()
            elif tool == "Soapfuse":
                pred_res_dict = self.get_soapfuse_results()
        # pylint exceptions:
        #        the caught exception is not further specified as several different exceptions may be raised during processing
        #        the type of exception is, however, unimportant for further processing, because any exception must be manually reviewed
        except Exception as ukn_err:  # pylint: disable=W0703
            logger.error("Couldn't fetch results from {0}, please check data in {1}. Error message: {2}".format(tool,
                                                                                                                output_folder_path,
                                                                                                                ukn_err))
            return (pred_res_dict, True)

        tool_res_file = os.path.join(output_folder_path, tool + "_res.csv")
        with open(tool_res_file, "w") as tool_outf:
            tool_outf.write("bpid;fusion_gene;breakpoint1;breakpoint2;junc_reads;span_reads;sample_id;tool\n")
            for key in pred_res_dict:
                tool_outf.write(key + ";" + ";".join(pred_res_dict[key]) + "\n")
        return (pred_res_dict, False)

    # pylint: enable=W0703

    # urla: naming is terribly complicated to understand and follow, changed every single var name!
    #       todo: the method seems to be very complicated for a relatively simple task => is there a more simple solution?
    # @param dict_of_fusion_results: This is a dictionary of dictionaries.
    #        For each fusion prediction tool, a dictionary is created with keys=fuid
    #        and value=fusion_info; these are themselfes organised in a dict
    #        with keys=fusion_tool and values=fusion_tool_dict
    def lookup_fusions_in_prediction(self, dict_of_fusion_dicts):
        """UKN"""
        # @param dict_of_fusionid_lists: This is a dictionary of lists.
        #        For each fusion tool (=keys of the dict) it contains a
        #        list (=value of the dict) the keys from the respective
        #        fusion tool dict
        dict_of_fusionid_lists = {}

        # for every tool in the results dict
        for fusion_tool in dict_of_fusion_dicts:
            # if no fusions were predicted by a tool, set an empty list
            # urla: is this actually required? it would anyway become an empty
            #       list in the for loop, wouldn't it?
            if len(dict_of_fusion_dicts[
                       fusion_tool]) == 0:  # this is a "number of elements" check => pylint: disable=C1801
                dict_of_fusionid_lists[fusion_tool] = []

            # for each fusion id (=keys of the respective fusion dict),
            # append them to a list in the new dict
            # urla: todo: the exception should only be raised, if the list is
            #             not existing (remove generalization)
            #       todo2: why not initializing an empty list for all in advance?
            #              i.e.: putting "dict_of_fusionid_lists[fusion_tool] = []"
            #              at the start of the for loop (would also eliminate the if len(...))
            for fusion_id in dict_of_fusion_dicts[fusion_tool]:
                try:
                    dict_of_fusionid_lists[fusion_tool].append(fusion_id)
                except KeyError:
                    print(
                        "Error when trying to append to list of tool {0}. Trying to create new list with {1} at start".format(
                            fusion_tool, fusion_id))
                    dict_of_fusionid_lists[fusion_tool] = [fusion_id]

        # create a list of unique fusion ids
        list_of_all_fusion_ids = []
        for tool in dict_of_fusionid_lists:
            list_of_all_fusion_ids += dict_of_fusionid_lists[tool]
        list_of_unique_fusion_ids = list(set(list_of_all_fusion_ids))

        # @param dict_of_found_uniq_fusions: A dictionary of lists, where keys
        #        are the unique fusion ids and values a list of booleans indicating
        #        whether or not the fusion was found in a fusion prediction tool
        dict_of_found_uniq_fusions = {}
        for uniq_fusion_id in list_of_unique_fusion_ids:
            # split the fusion id into [0]=gene1, [1]=breakpoint of gene1,
            # [2]=gene2, [3]=breakpoint of gene2
            uniq_fusion_id_split = uniq_fusion_id.split("_")

            list_of_found_fusion_booleans = []
            # for each fusion tool
            for fusion_tool in dict_of_fusion_dicts:
                found_fusion = False
                # for each fusion id in the list of fusions per tool
                for fusion_id in dict_of_fusionid_lists[fusion_tool]:
                    fusion_id_split = fusion_id.split("_")
                    # check if gene1 and gene2 of the unique fusion which is being tested
                    # is present in the current fusion (independent of the orientation)
                    # urla: this will lead to false positive results:
                    #       eg: fusion with gene1 "AB1" and gene2 "AB2" will match to fusion "AB11"-"AB22"
                    # if uniq_fusion_id_split[0] in fusion_id_split and uniq_fusion_id_split[2] in fusion_id_split:
                    # urla: possible solution:
                    # print(fusion_id_split)
                    # print(uniq_fusion_id_split)
                    if ((uniq_fusion_id_split[0] == fusion_id_split[0] and uniq_fusion_id_split[1] == fusion_id_split[
                        1]) or
                            (uniq_fusion_id_split[0] == fusion_id_split[1] and uniq_fusion_id_split[1] ==
                             fusion_id_split[0])):
                        found_fusion = True
                        break  # we don't need to look further if it was found at least once
                list_of_found_fusion_booleans.append(found_fusion)

            if sum(list_of_found_fusion_booleans) >= self.tool_num_cutoff:
                dict_of_found_uniq_fusions[uniq_fusion_id] = list_of_found_fusion_booleans
        return dict_of_found_uniq_fusions

    def run(self):
        """ asd """
        tool_state_path = os.path.join(self.fusion_output_path, "tool_state.csv")
        with open(tool_state_path, "w") as tool_state:
            tool_state.write("Sample ID")
            for tool in self.tools:
                tool_state.write(", {}".format(tool))
            tool_state.write("\n")

        detected_fusions_file = os.path.join(self.fusion_output_path, "Detected_Fusions.csv")
        with open(detected_fusions_file, "w") as fus_file:
            # write header
            fus_file.write("BPID;Fusion_Gene;Breakpoint1;Breakpoint2;Junction_Reads;Spanning_Reads;Sample;Tool\n")
            count_fusions = 0
            logger.debug("Generating Detected Fusions table")

            fusion_result_dict, results_with_errors_dict = self.concatenate_fusion_results(tool_state_path,
                                                                                           self.fusion_output_path)
            # print(len(fusion_result_dict))
            if sum(results_with_errors_dict.values()) == len(fusion_result_dict):
                logger.error("Fusion parsing failed completely. Revision required. Aborting.")
                sys.exit(1)
            elif sum(results_with_errors_dict.values()) != 0:
                logger.error("Results incomplete. Please make sure that all tools have run completely on dataset.")

            dict_of_boolean_list_of_found_uniq_fusions = self.lookup_fusions_in_prediction(
                fusion_result_dict)  # this is snake case, although too long... pylint: disable=C0103
            # for each unique fusion
            for uniq_fusion_id in dict_of_boolean_list_of_found_uniq_fusions:
                # for each fusion tool
                for fusion_tool_num_in_list, fusion_tool in enumerate(fusion_result_dict, 0):
                    # if the unique fusion was found in a fusion tool
                    if dict_of_boolean_list_of_found_uniq_fusions[uniq_fusion_id][fusion_tool_num_in_list]:
                        # for each fusion id of fusion tool X
                        for fusion_id in fusion_result_dict[fusion_tool]:
                            # if the fusion id of the tool matches the unique fusion id, write everything to file
                            # urla: it is probably better to iterate over all fusion, instead of using "get" on the dict
                            #       because a fusion gene can be called more than once with sligthly different breakpoints
                            if fusion_id == uniq_fusion_id:
                                count_fusions += 1
                                fus_file.write(
                                    uniq_fusion_id + ";" + ";".join(fusion_result_dict[fusion_tool][fusion_id]) + "\n")

        logger.info("Wrote {0} detected fusion genes to {1}.".format(count_fusions, detected_fusions_file))


def main():
    """Main"""
    parser = ArgumentParser(description='Extracts information on fusion genes')
    parser.add_argument('-i', '--input', dest='input', help='Specify the scratch folder of the sample.', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Specify the fusion output folder of the sample.',
                        required=True)
    parser.add_argument('-s', '--sample', dest='sample', help='Specify the sample to process.', required=True)
    parser.add_argument('-t', '--toolcutoff', dest='toolcutoff',
                        help='The minimal of tools that must support a fusion in order to be considered further.',
                        required=True)
    parser.add_argument('-f', '--fusionlist', dest='fusionlist',
                        help='A list of fusion prediction tools that where run on the sample.', required=True)
    parser.add_argument('-l', '--logger', dest='logger', help='Logging of processing steps', default="")
    args = parser.parse_args()

    res_parse_1 = FusionParser(args.input, args.output, args.sample, args.toolcutoff, args.fusionlist, args.logger)
    res_parse_1.run()


if __name__ == '__main__':
    main()
