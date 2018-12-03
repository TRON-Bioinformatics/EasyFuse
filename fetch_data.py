#!/usr/bin/env python

"""
Collect fusion prediction results,
annotate fusions,
run requantification and
combine all results

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""

from __future__ import print_function
import sys
import os
import time
import math

from argparse import ArgumentParser

from misc.config import Config
import misc.logger as Logger
import misc.io_methods as IOMethods

class ResultParser(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""
    def __init__(self, config, working_dir, prio_cutoff, stringency):
        """Parameter initiation and work folder creation."""
        self.cfg = cfg
        self.working_dir = working_dir
        self.prio_cutoff = prio_cutoff
        self.stringency = stringency
        self.fetchdata_dir = os.path.join(self.working_dir,"fetchdata_{} tool".format(self.prio_cutoff))
        self.sample_id = self.working_dir.rstrip("/").rstrip("/scratch").rsplit("/",1)[1]
#        self.sub_dict = {"orf":"ORF","GCN1L1":"GCN1","PTPLB":"HACD2"}
        self.chr_list = [
            "1", "2", "3", "4", "5",
            "6", "7", "8", "9", "10",
            "11", "12", "13", "14", "15",
            "16", "17", "18", "19", "20",
            "21", "22", "X", "Y"
        ]

        logfile = os.path.join(self.working_dir,"fusion.log")
        self.logger = logger.Logger(logfile)
        self.logger.setLevel(logger.DEBUG)
        self.logger.debug("Starting Result Parser")

        IOMethods.create_folder(self.fetchdata_dir)

    def get_fusioncatcher_results(self):
        summary_file = os.path.join(
            self.working_dir, 
            "fusion", 
            "fusioncatcher", 
            "summary_candidate_fusions.txt"
        )
        summary_file_detail = os.path.join(
            self.working_dir, 
            "fusion", 
            "fusioncatcher", 
            "final-list_candidate-fusion-genes.txt"
        )
        reciprocal_fusions = []
        with open(summary_file) as summary:
#            start = False
            for line in summary:
                if line.strip().startswith("*"):
                    fusion_gene = re.search("\*\s([\S]*)", line).group(2)
                    fp_write(fusion_gene + "\n")

                    if "reciprocal" in line:
                        reciprocal_fusions.append(fusion_gene.replace("--", "_").upper())
        fusion_map = {}
        with open(summary_file_detail) as summary_detail:
            # skip header line
            next(summary_detail)
            for line in summary_detail:
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

                fusion_gene = "_".join([elements[0], elements[1]]).upper()
                
                
                if fusion_gene in reciprocal_fusions:
                    fusion_gene = "_".join([elements[1], elements[0]]).upper()
                
                if elements[8].split(":")[0] not in self.chr_list or elements[9].split(":")[0] not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = "_".join([fg_split[0], elements[8], fg_split[1], elements[9]])

#                fusion_map[fgid] = [fusion_gene,up_gene_bp,dn_gene_bp,junc_reads_num,span_reads_num,self.sample_id,"fusioncatcher"]


                fusion_map[fgid] = [
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
    def get_starfusion_results(self):
        """Load and parse results from star-fusion"""
        summary_file = os.path.join(
            self.working_dir, 
            "fusion", 
            "starfusion", 
            "star-fusion.fusion_candidates.final"
        )
        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for line in summary:
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

                fusion_gene = elements[0].replace("--","_").upper()
                
                if elements[5].split(":")[0] not in self.chr_list or elements[7].split(":")[0] not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = "_".join([fg_split[0], up_gene_bp, fg_split[1], dn_gene_bp])

                fusion_map[fgid] = [
                    fusion_gene,  # fusion_gene
                    elements[5],  # up_gene_bp
                    elements[7],  # dn_gene_bp
                    elements[1],  # junc_reads_num
                    elements[2],  # span_reads_num
                    self.sample_id,
                    "Starfusion"
                ]
        return fusion_map


    def get_soapfuse_results(self):
        summary_file = ""
        folder_to_scan = os.path.join(
            self.working_dir,
            "fusion",
            "soapfuse",
            "final_fusion_genes"
        )
        for file in os.listdir(folder_to_scan):
            folder_path = os.path.join(folder_to_scan,file)
            if os.path.isdir(folder_path):
                for res in os.listdir(folder_path):
                    if res.endswith(".final.Fusion.specific.for.genes"):
                        summary_file = os.path.join(folder_path,res)
        if not summary_file:
            summary_file = os.path.join(
                self.working_dir,
                "fusion",
                "soapfuse",
                "final_fusion_genes",
                self.sample_id,
                self.sample_id + ".final.Fusion.specific.for.genes"
            )
        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for i,line in enumerate(summary):
                elements = line.rstrip().split("\t")
                # Currently relevant fields (marked *) in the output file are:
                # * elements[0] ~ Gene_1_symbol(5end_fusion_partner)
                # * elements[1] ~ chromosome of up stream fusion partner
                # * elements[2] ~ strand of up stream fusion partner
                # * elements[3] ~ genome junction position of up stream fusion partner
                #   elements[4] ~ location of up stream fusion partner's junction point
                # * elements[5] ~ Gene_2_symbol(3end_fusion_partner)
                # * elements[6] ~ chromosome of down stream fusion partner
                # * elements[7] ~ strand of down stream fusion partner
                # * elements[8] ~ genome junction position of down stream fusion partner
                # * elements[9] ~ location of down stream fusion partner's junction point
                # * elements[10] ~ number of span-reads (details in sample-ID.final.Span_reads
                # * elements[11] ~ number of junc-reads (details in sample-ID.final.Junc_reads)
                #   elements[12] ~ classification of fusions
                #   elements[13] ~ whether down stream fusion partner is frame-shift or in-frame-shift

                fusion_gene = "_".join([elements[0], elements[5]]).upper()

                up_gene_bp = ":".join([elements[1], elements[3], elements[2]])
                dn_gene_bp = ":".join([elements[6], elements[8], elements[7]])
                
                if elements[1] not in self.chr_list or elements[6] not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = "_".join([fg_split[0], up_gene_bp, fg_split[1], dn_gene_bp])
                fusion_map[fgid] = [
                    fusion_gene,
                    up_gene_bp,
                    dn_gene_bp,
                    elements[11],
                    elements[10],
                    self.sample_id,
                    "soapfuse"
                ]
        return fusion_map

    def get_mapsplice_results(self):
        """Load and parse results from mapsplice"""
        summary_file = os.path.join(
            self.working_dir,
            "fusion",
            "mapsplice",
            "fusions_well_annotated.txt"
        )

        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for line in summary:
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
                
                
                gene_name_1 = elements[60].split(",")[0]
                gene_name_2 = elements[61].split(",")[0]
                fusion_gene = "_".join([gene_name_1, gene_name_2]).upper()

                chrs = elements[0].split("~")
                up_gene_bp = ":".join([chrs[0], elements[1], elements[5][0]])
                dn_gene_bp = ":".join([chrs[1], elements[2], elements[5][1]])
                
                if chrs[0] not in self.chr_list or chrs[1] not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = "_".join([fg_split[0], up_gene_bp, fg_split[1], dn_gene_bp])

                fusion_map[fgid] = [
                    fusion_gene,
                    up_gene_bp,
                    dn_gene_bp,
                    elements[4],
                    elements[27],     # previous elements[19]
                    self.sample_id,
                    "mapsplice"
                ]
        return fusion_map


    # starchip - results file is "starchip.summary"
    def get_starchip_results(self):
        """Load and parse results from starchip"""
        starchip_predict = os.path.join(
            self.scratch_path, 
            "fusion", 
            "starchip", 
            "starchip.summary"
        )
        fusion_map = {}
        with open(starchip_predict, "r") as prediction:
#            with open(os.path.join(self.fusion_output_path, "Starchip_res_toFusionInspector.txt"), "w") as fp_write:
            next(prediction) # skip header line
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
                fusion_gene = "_".join([elements[5], elements[7]]).upper()
#                fp_write.write("{0}\n".format(fusion_gene.replace("_", "--")))
                    
                # check whether fusion gene is not on primary chr
                if elements[0].split(":")[0] not in self.chr_list or elements[1].split(":")[0] not in self.chr_list:
                    continue
                fgid = "_".join([elements[5], elements[0], elements[7], elements[1]])
                        
                fusion_map[fgid] = [
                    fusion_gene,  # fusion_gene
                    elements[0],  # up_gene_bp
                    elements[1],  # dn_gene_bp
                    elements[3],  # junc_reads_num
                    elements[2],  # span_reads_num
                    self.sample_id,
                    "Starchip"
                ]
        return fusion_map

    def get_jaffa_results(self):
        summary_file = os.path.join(
            self.working_dir, 
            "fusion", 
            "jaffa", 
            "jaffa_results.csv"
        )

        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for line in summary:
                elements = line.rstrip().split(",")
                elements = [ele.strip("\"") for ele in elements]
                fusion_gene = elements[1].replace(":","_").upper()

                up_gene_bp = ":".join([elements[2], elements[3], elements[4]])
                dn_gene_bp = ":".join([elements[5], elements[6], elements[7]])
                
                if elements[2] not in self.chr_list or elements[5] not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = "_".join([fg_split[0], up_gene_bp, fg_split[1], dn_gene_bp])

                fusion_map[fgid] = [
                    fusion_gene,
                    up_gene_bp,
                    dn_gene_bp,
                    elements[9],
                    elements[10],
                    self.sample_id,
                    "jaffa"
                ]
        return fusion_map

    def get_infusion_results(self):
        summary_file = os.path.join(self.working_dir,"fusion","infusion","fusions.detailed.txt")

        fusion_map = {}
        with open(summary_file) as summary:
            next(summary)
            for line in summary:
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
                
#                fusion_gene = (elements[21] + "_" + elements[27]).upper().replace(";","-")
                gene_name_1 = self.get_fusion_gene_name_infusion(elements[21], elements[22])
                gene_name_2 = self.get_fusion_gene_name_infusion(elements[27], elements[28])
                fusion_gene = "_".join([gene_name_1, gene_name_2])

                # To-Do: correct none-type strand ("." -> "+"/"-")
                up_gene_bp = ":".join([elements[1], elements[2], elements[23]])
                dn_gene_bp = ":".join([elements[4], elements[5], elements[29]])
                
                if elements[1] not in self.chr_list or elements[4] not in self.chr_list:
                    continue

                fg_split = fusion_gene.split("_")
                fgid = "_".join([fg_split[0], up_gene_bp, fg_split[1], dn_gene_bp])

                fusion_map[fgid] = [
                    fusion_gene,
                    up_gene_bp,
                    dn_gene_bp,
                    elements[7],
                    elements[8],
                    self.sample_id,
                    "infusion"
                ]
        return fusion_map

    def get_fusion_gene_name_infusion(gene_id, transcript_list_field):
        """Helper method for infusion data parsing. Returns the most frequently listed fusion id in a list of fusion ids"""
        # if only 1 gene id is returned as fusion partner, return this
        if not ";" in gene_id:
            return gene_id.replace("_","-")

        # if 2 or more gene ids are given, take the most frequent from the transcript list
        gid_dict = {}
        for gid in transcript_list_field.split(";"):
            id_test = gid.split(":")[0][:-4]
            if not id_test in gid_dict:
                gid_dict[id_test] = 0
            gid_dict[id_test] += 1
                
        best_hit = max(gid_dict.iteritems(), key=operator.itemgetter(1))[0]
                
        return best_hit.replace("_","-")
                

    # pizzly - results file is "kallizzy.json.txt"
    def get_pizzly_results(self):
        """Load and parse results from pizzly"""
        pizzly_predict = os.path.join(self.scratch_path, "fusion", "pizzly", "kallizzy.json.txt")
        fusion_map = {}
        with open(pizzly_predict, "r") as prediction:
#            with open(os.path.join(self.fusion_output_path, "Pizzly_res_toFusionInspector.txt"), "w") as fp_write:
            next(prediction) # skip header line
            for line in prediction:
                elements = line.rstrip().split("\t")
                # Currently relevant fields (marked *) in the output file are:
                # * elements[0] ~ geneA.name
                #   elements[1] ~ geneA.id
                # * elements[2] ~ geneB.name
                #   elements[3] ~ geneB.id
                # * elements[4] ~ paircount
                # * elements[5] ~ splitcount
                #   elements[6] ~ transcripts.list
                
                # Pizzly is a overpredictor with a high FP ratio in its current version
                # Therefore, only events supported by at least (paircount + splitcount >=3) are considered
                if int(elements[4]) + int(elements[5]) < 3:
                    continue
                fusion_gene = "{0}_{1}".format(elements[0], elements[2]).upper()
                # fp_write.write("{0}--{1}\n".format(elements[0], elements[2]))

                # check whether fusion gene is not on primary chr - not possible for pizzly as exact breakpoint cannot directly be determined from pizzly output
                # if elements[0].split(":")[0] not in self.chr_list or elements[1].split(":")[0] not in self.chr_list:
                # continue
                fgid = elements[0] + "_1:100:+_" + elements[2] + "_2:100:+"
                        
                fusion_map[fgid] = [
                    fusion_gene,  # fusion_gene
                    "1:100:+",  # up_gene_bp
                    "2:100:+",  # dn_gene_bp
                    elements[5],  # junc_reads_num
                    elements[4],  # span_reads_num
                    self.sample_id,
                    "Pizzly"
                ]
        return fusion_map


    def generate_matrix(self, res):
        keys = {}
        for tool in res:
            if len(res[tool]) == 0:
                keys[tool] = []

            for key in res[tool]:
                try:
                    keys[tool].append(key)
                except:
                    keys[tool] = [key]

        total_fusion_genes = []
        for tool in keys:
            total_fusion_genes += keys[tool]
        fusion_genes = list(set(total_fusion_genes))

        matrix = {}
        for key in fusion_genes:
            key_split = key.split("_")

            t = []
            for tool in res:
                found = False
                for ele in keys[tool]:
                    ele_split = ele.split("_")
                    if key_split[0] in ele_split and key_split[2] in ele_split:
                        found = True

                t.append(found)

            sum_tools = sum(t)
            if sum_tools >= self.prio_cutoff:
                matrix[key] = t
        return matrix

    def concatenate_fusion_results(self, tool_state_path):
        """Return tuple of (dict of results dicts, dict of errors) and writes error/pass to summary file"""
        fusion_res_dict = {}
        results_with_errors_dict = {}
        tools = self.cfg.get('general','tools').split(",")

        self.logger.info("Processing " + self.sample_id)
        output_folder_path = os.path.join(self.fetchdata_dir, "tool_res")
        IOMethods.create_folder(output_folder_path)
        sample_string = self.sample_id
        for tool in tools:
            fusion_res_dict[tool], results_with_errors_dict[tool] = self.get_tool_results(output_folder_path, tool)

            if results_with_errors_dict[tool]:
                sample_string += ";0"
            else:
                sample_string += ";1"
        sample_string.rstrip(";")

        with open(tool_state_path, "a") as outf:
            outf.write(sample_string + "\n")

        return (fusion_res_dict, results_with_errors_dict)

    def get_tool_results(self, output_folder_path, tool):
        """Return tuple of (dict of results from individual fusion tools, error type)"""
        pred_res_dict = {}
        self.logger.info("Parsing results for " + tool)
        try:
            if tool == "fusioncatcher":
                pred_res_dict = self.get_fusioncatcher_results()
            elif tool == "starfusion":
                pred_res_dict = self.get_starfusion_results()
            elif tool == "soapfuse":
                pred_res_dict = self.get_soapfuse_results()
            elif tool == "mapsplice":
                pred_res_dict = self.get_mapsplice_results()
            elif tool == "jaffa":
                pred_res_dict = self.get_jaffa_results()
            elif tool == "infusion":
                pred_res_dict = self.get_infusion_results()
        except Exception as ukn_err:
            self.logger.error("Couldn't fetch results from {0}, please check data in {1}. Error message: {2}".format(tool, output_folder_path, ukn_err))
            return (pred_res_dict, True)

        tool_res_file = os.path.join(output_folder_path, tool + "_res.csv")
        with open(tool_res_file, "w") as tool_outf:
            tool_outf.write("fgid;fusion_gene;breakpoint1;breakpoint2;junc_reads;span_reads;sample_id;tool\n")
            for key in pred_res_dict:
                tool_outf.write("{};{}\n".format(key, ";".join(res[key])))

        return (pred_res_dict, False)

    def run(self):
        tools = self.cfg.get('general', 'tools').split(",")
        tool_state_path = os.path.join(self.fetchdata_dir, "tool_state.csv")
        with open(tool_state_path, "w") as tool_state:
            tool_state.write("Sample ID")
            for tool in tools:
                tool_state.write(", {}".format(tool))
            tool_state.write("\n")


        res, err = self.concatenate_fusion_results(tool_state_path)

        detected_file = os.path.join(self.fetchdata_dir, "Detected_Fusions.csv")
        outf = open(detected_file, "w")

        outf.write("FGID;Fusion_Gene;Breakpoint1;Breakpoint2;Junction_Reads;Spanning_Reads;Sample;Tool\n")
        c = 0
        self.logger.debug("Generating Detected Fusions table")

        if sum(err.values()) != 0:
            sys.exit(1)

        matrix = self.generate_matrix(res)
        for key in matrix:
            for i, tool in enumerate(res):
                if matrix[key][i]:
                    for f_key in res[tool]:
                        if f_key == key:
                            c += 1
                            outf.write("{};{}\n".format(key, ";".join(res[tool][f_key])))

        outf.close()

        self.logger.debug("Wrote {} lines.".format(c))
        self.logger.info("Detected fusion genes table created.")

        if incomplete and not self.stringency:
            self.logger.error("Results incomplete. Aborting. Please make sure that all tools have run completely on dataset.")
            sys.exit(1)

        context_seqs = os.path.join(self.fetchdata_dir, "Context_seqs.csv")

#    cmd = "%s {} {} {}" % (cfg.get('commands','fetch_context_cmd'),detected_file,fetchdata_dir,organism)
        cmd = "{} -i {} -o {} -f {} -e {} --context_seq_len {}".format(self.cfg.get('commands','fetch_context_cmd'), detected_file, context_seqs, self.cfg.get('bowtie_indexes','bowtie_index_hg38_fasta'), self.cfg.get('references','ensembl_genes_hg38_csv'), self.cfg.get('general', 'context_seq_len'))

        if os.path.exists(detected_file) and not os.path.exists(context_seqs):

            self.logger.info("Generating context sequences out of detected fusions.")
            self.logger.debug("Command was: " + cmd)
            proc = subprocess.Popen(["/bin/bash", "-c", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (stdoutdata,stderrdata) = proc.communicate()
            if proc.returncode != 0:
                self.logger.error(stderrdata)

#        seq2hla_classI_res = os.path.join(self.working_dir,"seq2hla","run-ClassI-class.HLAgenotype4digits")
#        seq2hla_classII_res = os.path.join(self.working_dir,"seq2hla","run-ClassII.HLAgenotype4digits")

#        cmd = "{} {} {} {} {}" % (self.cfg.get('commands','mhc_prediction_cmd'),context_seqs,seq2hla_classI_res,seq2hla_classII_res,self.fetchdata_dir)

#        classI_pred = os.path.join(self.fetchdata_dir,"hla_classI.csv")
#        classII_pred = os.path.join(self.fetchdata_dir,"hla_classII.csv")

#        if not os.path.exists(classI_pred) and not os.path.exists(classII_pred) and os.path.exists(context_seqs):
#            self.logger.info("Predicting MHC eptitopes for peptides.")
#            self.logger.debug("Command was: " + cmd)
#            proc = subprocess.Popen(["/bin/bash","-c",cmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#            (stdoutdata,stderrdata) = proc.communicate()
#            if proc.returncode != 0:
#                self.logger.error(stderrdata)

        file_list = []

        file_r1 = os.path.join(self.working_dir,"skewer","out_file-trimmed-pair1.fastq.gz")
        file_r2 = os.path.join(self.working_dir,"skewer","out_file-trimmed-pair2.fastq.gz")
        if not os.path.exists(file_r1) and not os.path.exists(file_r2):
            file_r1 = os.path.join(self.working_dir,"out_file-trimmed-pair1.fastq.gz")
            file_r2 = os.path.join(self.working_dir,"out_file-trimmed-pair2.fastq.gz")
        file_list.append(file_r1)
        file_list.append(file_r2)


#    requantification_dir = os.path.join(working_dir,"requantification")
        cmd = "{} -i {} -f {} -o {} -p prod".format(self.cfg.get('commands','custom_trans_cmd')," ".join(file_list),context_seqs,self.fetchdata_dir)
    
        classification_file = os.path.join(self.fetchdata_dir,"Classification.csv")
#    classification_file_extend = os.path.join(working_dir,"Classification_extend.csv")

        if not os.path.exists(classification_file) and os.path.exists(context_seqs):
            self.logger.info("Creating custom transcriptome and quantifying samples.")
            self.logger.debug("Command was: " + cmd)
            proc = subprocess.Popen(["/bin/bash","-c",cmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            (stdoutdata,stderrdata) = proc.communicate()
            if proc.returncode != 0:
                self.logger.error(stderrdata)

        if os.path.exists(context_seqs):
            self.logger.info("Adding fusion expression to table")
            self.add_fusion_expression(context_seqs,quant_file)

####
#    res_parse.add_gene_expression(args.input,args.output,quant_file)

#    print "Adding fusion frequency to table"
#    res_parse.add_frequency(args.input,outfile+".state.expression.fusion")
    
#    print "Adding exon boundaries to table"
#    res_parse.add_exon_boundaries(context_seqs,outfile+".state.expression.fusion.freq")

#    print "Filtering final table"
#    res_parse.filter_final_table(outfile+".state.expression.fusion.freq.exbnd.csv")
    
#    primerblat_results = os.path.join(args.output,"fusion_primerblat.csv")
#    primerblat_err = os.path.join(args.output,"fusion_primerblat.err")

#    cmd = "python /kitty/code/primerblat/primerblat.py --primer3_target=80,40 --primer3_size=100-150  --column_values=1 --inputdelim=";"  --debug_p3 --delim=";" --dec="," {} /kitty/data/human/hg38.2bit > {} 2> {}" % (context_seqs,primerblat_results,primerblat_err)
#    log.info("Generating primers out of context sequences.")
#    log.debug("Command was: " + cmd)
#    proc = subprocess.Popen(["/bin/bash","-c",cmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#    (stdoutdata,stderrdata) = proc.communicate()
#    if proc.returncode != 0:
#        log.error(stderrdata)        
                

def main():
    parser = ArgumentParser(description='Extracts information on fusion genes')
    parser.add_argument('-i', '--input', dest='input', help='Specify sample working directory.',required=True)
    parser.add_argument('-c', '--config', dest='config', help='Specify config file.',default="")
#    parser.add_argument('-g', '--gene-symbols', dest='gene_symbols', help='Specify hugo to enembl gene symbol mapping file',default='/kitty/data/human/ensembl/GRCh38.86/Homo_sapiens.GRCh38.86.HGNC2ENS.csv')
    parser.add_argument('-s', '--species', dest='species', choices=['human','mouse'], help='Specify species to be processed',default='human')
    parser.add_argument('-v', '--stringent', dest='stringency', action='store_true', help='Special case where not all tools have results')

    args = parser.parse_args()

    if args.config == "":
        cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)),'config.ini'))
    else:
        cfg = Config(args.config)
    working_dir = args.input
    organism = args.species

    if not os.path.exists(working_dir) or not os.path.isdir(working_dir):
        sys.exit(1)
    
    script_call = "python " + os.path.realpath(__file__) + " " + " ".join(sys.argv[1:])

    outf = open(os.path.join(working_dir,"fetch_data.sh"),"w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

    
    print "Running with 1 tool cutoff"
    try:
        res_parse_1 = ResultParser(cfg,working_dir,1,args.stringency)
        res_parse_1.run()
    except:
        pass

    print "Running with 2 tool cutoff"
    try:
        res_parse_2 = ResultParser(cfg,working_dir,2,args.stringency)
        res_parse_2.run()
    except:
        pass

            
                
if __name__ == '__main__':
    main()
