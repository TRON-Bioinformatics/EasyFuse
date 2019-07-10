#!/usr/bin/env python

"""
Annotate predicted fusion genes. This involves:
    1) Identify transcripts possibly involved in the fusion
    2) Determine the translation frame at the breakpoint position
    3) Reconstruct possible fusion transcripts as well as the respective wildtype background
    4) Translate transcripts
    5) Select regions of interest for neopeptide searches
    6) Collect information on the kind of fusion which might be relevant for downstream selection processes
    6) Filter likely artifacts based on the quality of the annotated transcripts
@author: BNT (URLA)
@version: 20190704
"""
from __future__ import print_function, division
import argparse
import gffutils
import xxhash
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class FusionAnnotation(object):
    """ Annotation of predicted fusion genes soley based on the breakpoint information """
    def __init__(self, db_file, bp_file, tsl_file):
        """Parameter initialization"""
        self.db = self.load_featuredb(db_file)
        self.bp_dict = self.breakpoints_to_dict(bp_file)
        self.tsl_dict = self.load_tsl_data(tsl_file)
        self.cds_seq_dict = {}
        self.trans_flags_dict = {}

    def load_featuredb(self, db_in_file):
        """ Perform a simple check to confirm that the database input is a database """
        # urla - note: I would like to have a better check implemented before trying to load the db but currently don't know how this should look like
        if not db_in_file.endswith("db"):
            print("Not a valid database file: {}".format(db_in_file))
        # return the db
        return gffutils.FeatureDB(db_in_file)

    def breakpoints_to_dict(self, bp_in_file):
        """ Read fgid and breakpoints from Detected_Fusions.csv and store the info in a dict """
        bp_dict = {}
        with open(bp_in_file, "r") as bpin:
            next(bpin) # skip header
            for line in bpin:
                # for the annotaion, we just need the fgid and the breakpoints
                fgid, _, bp1, bp2, _, _, _, _ = line.split(";")
                if not fgid in bp_dict:
                    bp_dict[fgid] = (bp1, bp2)
        return bp_dict

    def load_tsl_data(self, tsl_in_file):
        """ Load data created by gtf2tsl.py into a dict """
        tsl_dict = {}
        with open(tsl_in_file, "r") as tslin:
            next(tslin) # skip header
            for line in tslin:
                # for the annotation, we just need the transcript id and the tsl
                trans, _, _, tsl = line.strip().split("\t")
                if not trans in tsl_dict:
                    tsl_dict[trans] = tsl
        return tsl_dict

    def get_tsl(self, trans_id):
        """ Get the transcript support level (tsl) from the pre-loaded dict """
        # urla - note: the tsl is not available for grch37 data (at least not v87/90) and I'm therefore using the tsl info from grch38 v90
        #              With the newer data, however, it is possible that some transcript ids are not in the database (because the respective
        #              transcript annotations were deleted or modified significantly).
        # urla - note: Although the "correct" tsl for a missing trans_id would be "NA", I'm using "6" in order to allow differentiation from
        #              from existing, "normal" tsl's which range from 1-5 and are NA for transcripts of unknown support
        if not trans_id in self.tsl_dict:
            return "6"
        return self.tsl_dict[trans_id]

    def fgid2fusiongene(self, fgid):
        """ Create the fusion gene string from the fgid """
        fg1, _, fg2, _ = fgid.split("_")
        return "_".join([fg1, fg2])

    def define_type(self, cis_near_distance, bp1_chr, bp1_pos, bp1_strand, bp2_chr, bp2_pos, bp2_strand):
        """ Define the fusion type based on the location and orientation of the fusion partners """
        # Fusion type:
        # -	"trans" (bp of fusion partners on different chromosomes, but with same strand)
        # -	"trans_inv" (bp of fusion partners on different chromosomes and with different strand)
        # -	"cis_inv" (bp of fusion partners on same chromosome, but with different strand)
        # -	"cis_trans" (bp of fusion partners on same chromosome and with same strand, but 2nd fusion partner is before (strandwise) 1st fusion partner)
        # -	"cis_far" (bp of fusion partners on same chromosome, same strand, 1st before 2nd, but more than 1Mb apart from each other (genomic distance))
        # -	"cis_near" (bp of fusion partners on same chromosome, same strand, 1st before 2nd and closer than 1Mb apart from each other (genomic distance))
        if bp1_chr == bp2_chr:
            if bp1_strand == bp2_strand:
                if bp1_strand == "+":
                    if bp2_pos - bp1_pos < 0:
                        return "cis_trans"
                    if bp2_pos - bp1_pos > cis_near_distance:
                        return "cis_far"
                    return "cis_near"
                else:
                    if bp1_pos - bp2_pos < 0:
                        return "cis_trans"
                    if bp1_pos - bp2_pos > cis_near_distance:
                        return "cis_far"
                    return "cis_near"
            return "cis_inv"
        else:
            if bp1_strand == bp2_strand:
                return "trans"
            return "trans_inv"

    def get_bp_overlapping_features(self, bp_chr, bp_pos):
        """ Get two lists with breakpoint overlapping exons and cds from the database """
        bp_exons = []
        bp_cds = []
        for efeature in self.db.region(region="{}:{}-{}".format(bp_chr, bp_pos - 1, bp_pos + 1), completely_within=False, featuretype = ("exon", "CDS")):
            if efeature.featuretype == "exon":
                bp_exons.append(efeature)
            elif efeature.featuretype == "CDS":
                bp_cds.append(efeature)
        # correct positions in such a way that the same position in the exon/CDS list represents the same parental transcript
        # exons are the scaffold as there can be exons w/o cds but not vv
        exon_transcripts = [x.attributes["Parent"][0] for x in bp_exons]
        cds_transcripts = [x.attributes["Parent"][0] for x in bp_cds]
        bp_cds = [bp_cds[cds_transcripts.index(x)] if x in cds_transcripts else "" for x in exon_transcripts]
        return bp_exons, bp_cds

    def get_boundary(self, bp_pos, efeature):
        """ Check whether the breakpoint position is on an exon or CDS boundary """
        if efeature == "":
            return "NA"
        if bp_pos == efeature.start:
            return "left_boundary"
        if bp_pos == efeature.stop:
            return "right_boundary"
        if efeature.start < bp_pos < efeature.stop:
            return "within"
        return "outside"

    def get_boundary2(self, boundary1, boundary2, strand1, strand2):
        """ Check whether boundaries of both fusion partners are in the correct orientation """
        # urla - note: "both" is returned to be compatible with the existing pipeline although the name is not very meaningful here
        if strand1 == "+":
            if strand2 == "+":
                if boundary1 == "right_boundary" and boundary2 == "left_boundary":
                    return "both"
            else:
                if boundary1 == "right_boundary" and boundary2 == "right_boundary":
                    return "both"
        else:
            if strand2 == "+":
                if boundary1 == "left_boundary" and boundary2 == "left_boundary":
                    return "both"
            else:
                if boundary1 == "left_boundary" and boundary2 == "right_boundary":
                    return "both"
        return "no_match"

    def get_parents(self, efeature):
        """ Get transcript and gene parents from exon feature """
        trans_id = ""
        trans_biotype = ""
        gene_id = ""
        gene_name = ""
        gene_biotype = ""
        description = ""
        for count_parents, parent in enumerate(self.db.parents(efeature.id), 1):
            if parent.id.startswith("transcript:"):
                trans_id = parent.id
                trans_biotype = parent.attributes["biotype"][0]
            elif parent.id.startswith("gene:"):
                gene_id = parent.id
                gene_name = parent.attributes["Name"][0]
                gene_biotype = parent.attributes["biotype"][0]
                try:
                    description = parent.attributes["description"][0].replace(";"," ")
                except KeyError:
                    description = "NA"
        if count_parents > 2:
            print("Warning: {0} has more than two parents. Please check that gene_id {1} and trans_id {2} are correct for the exon at position {3}-{4}.".format(
                    efeature.id, gene_id, trans_id, efeature.start, efeature.stop))
        return trans_id, trans_biotype, gene_id, gene_name, gene_biotype, description

    def get_wt_codings(self, trans_id):
        """ Get matched exons and cds region per transcript id """
        exon_pos_list = []
        cds_pos_list = []
        cds_frame_list = []
        for feature in self.db.children(trans_id, featuretype = ("exon", "CDS"), order_by = "start", reverse = False):
            if feature.featuretype == "exon":
                exon_pos_list.append((feature.start, feature.stop))
            elif feature.featuretype == "CDS":
                cds_pos_list.append((feature.start, feature.stop))
                cds_frame_list.append(feature.frame)
        # correct list positions in such a way that the a cds is enclosed by an exon or not all available
        # exons are the scaffold as there can be exons w/o cds but not vv
        cds_pos_list2 = []
        cds_frame_list2 = []
        exon_has_cds = False
        for exon_start, exon_stop in exon_pos_list:
            exon_has_cds = False
            for cds_counter, (cds_start, cds_stop) in enumerate(cds_pos_list, 0):
                if cds_start >= exon_start and cds_stop <= exon_stop:
                    cds_pos_list2.append((cds_start, cds_stop))
                    if cds_start > exon_start and cds_stop < exon_stop:
                        self.trans_flags_dict[trans_id].add("cds != exon at no {}".format(cds_counter))
                    cds_frame_list2.append(cds_frame_list[cds_counter])
                    exon_has_cds = True
            if not exon_has_cds:
                cds_pos_list2.append("NA")
                cds_frame_list2.append("NA")
        return exon_pos_list, cds_pos_list2, cds_frame_list2

    def check_exon_overlap(self, wt1_exon_pos_list, wt2_exon_pos_list, fusion_type):
        """ Returns true if exons overlap and false otherwise """
        # exon position lists are always sorted and we therefore don't need any information on the strand
        # gene on different chromosomes
        if fusion_type.startswith("trans"):
            return False
        # end of wt1 locus must be < start of wt2 locus or start of wt1 locus > end of wt2 locus
        if wt1_exon_pos_list[-1][1] < wt2_exon_pos_list[0][0] or wt1_exon_pos_list[0][0] > wt2_exon_pos_list[-1][1]:
            return False
        return True

    def get_frame(self, bp_pos, cds_pos_list, cds_frame_list, strand):
        """ Get/Calculate the frame at the start of the cds and at the breakpoint """
        # urla - note: the frame annotation is imho a little confusing in ensembl and defined as follows:
        #              0: The next full codon (i.e. 3bp codon) starts 0 bases from the current position => this is the first (0) base of a codon
        #              1: The next full codon (i.e. 3bp codon) starts 1 base from the current position => this is the third (2) base of a codon
        #              2: The next full codon (i.e. 3bp codon) starts 2 bases from the current position => this is the second (1) base of a codon
        pos_frame_idx = [0, 2, 1, 0, 2, 1]
#        neg_frame_idx = [2, 0, 1, 2, 0, 1]

        # if the frame at bp is non-determinable, it is set to -1
        frame_at_start = -1
        frame_at_bp = -1
        # zip pos and frame while excluding NAs
        pos_frame_list = [(x,int(y)) for x,y in zip(cds_pos_list, cds_frame_list) if not x == "NA"]
        # find breakpoint cds and get the frame at the breakpoint
        for cds_pos, cds_frame in pos_frame_list:

#            if cds_pos[0] <= bp_pos <= cds_pos[1]:
#                tmp_frame_idx = pos_frame_idx[pos_frame_idx.index(cds_frame):pos_frame_idx.index(cds_frame)+3]
#                frame_at_bp = tmp_frame_idx[(bp_pos - (cds_pos[0] - 1)) % 3]
#                break

            tmp_frame_idx = pos_frame_idx[pos_frame_idx.index(cds_frame):pos_frame_idx.index(cds_frame)+3]
#            if (strand == "+" and bp_pos == cds_pos[0]) or (strand == "-" and bp_pos == cds_pos[1]):
#                frame_at_bp = cds_frame
#            elif cds_pos[0] < bp_pos <= cds_pos[1]:
#                frame_at_bp = tmp_frame_idx[(bp_pos - (cds_pos[0] - 1)) % 3]

            if strand == "+":
#                tmp_frame_idx = pos_frame_idx[pos_frame_idx.index(cds_frame):pos_frame_idx.index(cds_frame)+3]
#                frame_at_bp = tmp_frame_idx[(bp_pos - (cds_pos[0] - 1)) % 3]
                if bp_pos == cds_pos[0]:
                    frame_at_bp = cds_frame
#                    break
                elif cds_pos[0] < bp_pos <= cds_pos[1]:
                    frame_at_bp = tmp_frame_idx[(bp_pos - (cds_pos[0] - 1)) % 3]
#                    break
#                    frame_at_bp = pos_frame_idx[((bp_pos + 1) - cds_pos[0]) % 3 + cds_frame]
            else:
#                tmp_frame_idx = neg_frame_idx[pos_frame_idx.index(cds_frame):pos_frame_idx.index(cds_frame)+3]

                if bp_pos == cds_pos[1]:
                    frame_at_bp = cds_frame
#                    break
                elif cds_pos[0] <= bp_pos < cds_pos[1]:
                    frame_at_bp = tmp_frame_idx[(cds_pos[1] - (bp_pos - 1)) % 3]
#                    break
#                    frame_at_bp = neg_frame_idx[(cds_pos[1] - (bp_pos + 1)) % 3 + cds_frame]
        # get the starting frame of the cds
        if not len(pos_frame_list) == 0:
            if strand == "+":
                frame_at_start = pos_frame_list[0][1]
            else:
#                frame_at_start = frame_idx[(pos_frame_list[-1][0][1] - pos_frame_list[-1][0][0]) % 3 + pos_frame_list[-1][1]]
                frame_at_start = pos_frame_list[-1][1]
        # return frame at the beginning of the first cds and at the breakpoint position, both corrected according the the reading strand
        return frame_at_start, frame_at_bp

    def get_frame2(self, frame1, frame2):
        """ Combine frame info from bp1 and 2 """
        if frame1 == -1:
            return "no_frame"
        if frame2 == -1:
            return "neo_frame"
        if frame1 == frame2:
            return "in_frame"
        return "out_frame"

    def get_fusion_feature_list(self, bp1_pos, bp2_pos, bp1_strand, bp2_strand, bp1_feature_pos_list, bp2_feature_pos_list):
        """ Based on the breakpoints, the strand and the complete feature positions of the involved genes, return only those feature positions which will remain in the fusion. """
        bp1_feature_fus_list = []
        bp2_feature_fus_list = []

        # get fusion partner 1 cds
        if not len(bp1_feature_pos_list) == 0:
            if bp1_strand == "+": # on the "+" strand, we need everything LEFT of the bp for fusion gene partner 1
                bp1_feature_fus_list = [feature for feature in bp1_feature_pos_list if bp1_pos >= feature[0]]
                if not len(bp1_feature_fus_list) == 0:
                    if bp1_feature_fus_list[-1][0] == bp1_pos: # delete the last cds (i.e. last before the bp) if it starts with the breakpoint
                        bp1_feature_fus_list = bp1_feature_fus_list[:-1]
                    else: # or shorten the cds to the bp if it is within it
                        bp1_feature_fus_list[-1] = (bp1_feature_fus_list[-1][0], bp1_pos)
            else: # on the "-" strand, we need everything RIGHT of the bp for fusion gene partner 1
                bp1_feature_fus_list = [feature for feature in bp1_feature_pos_list if bp1_pos <= feature[1]]
                if not len(bp1_feature_fus_list) == 0:
                    if bp1_feature_fus_list[0][1] == bp1_pos: # delete first cds (i.e. first after the bp) if it starts with the breakpoint
                        bp1_feature_fus_list = bp1_feature_fus_list[1:]
                    else:
                        bp1_feature_fus_list[0] = (bp1_pos, bp1_feature_fus_list[0][1])

        # get fusion partner 2 cds
        if not len(bp2_feature_pos_list) == 0:
            if bp2_strand == "+": # on the "+" strand, we need everything RIGHT of the bp for fusion gene partner 2
                bp2_feature_fus_list = [feature for feature in bp2_feature_pos_list if bp2_pos <= feature[1]]
                if not len(bp2_feature_fus_list) == 0:
                    if bp2_feature_fus_list[0][1] == bp2_pos: # delete first cds (i.e. first after the bp) if it starts with the breakpoint
                        bp2_feature_fus_list = bp2_feature_fus_list[1:]
                    else:
                        bp2_feature_fus_list[0] = (bp2_pos, bp2_feature_fus_list[0][1])
            else: # on the "-" strand, we need everything LEFT of the bp for fusion gene partner 2
                bp2_feature_fus_list = [feature for feature in bp2_feature_pos_list if bp2_pos >= feature[0]]
                if not len(bp2_feature_fus_list) == 0:
                    if bp2_feature_fus_list[-1][0] == bp2_pos: # delete the last cds (i.e. last before the bp) if it starts with the breakpoint
                        bp2_feature_fus_list = bp2_feature_fus_list[:-1]
                    else: # or shorten the cds to the bp if it is within it
                        bp2_feature_fus_list[-1] = (bp2_feature_fus_list[-1][0], bp2_pos)
        return bp1_feature_fus_list, bp2_feature_fus_list

    def fill_seq_lookup_dict(self, chrom, cds_list):
        """ Create/fill a dictionary of genomic position tuple with chromosomes as keys and a list of two lists as values.
            Genomic coordinates are loaded into the first list and the second one is later complemented by the respective sequence"""
        if not chrom in self.cds_seq_dict:
            self.cds_seq_dict[chrom] = [[], []]
        for cds_coord in cds_list:
            if not cds_coord in self.cds_seq_dict[chrom][0]:
                self.cds_seq_dict[chrom][0].append(cds_coord)

    def get_stranded_seq(self, sequence, strand):
        """ Return the reverse complement of a sequence if the strand is negative and the unchanged sequence otherwise """
        if strand == "-":
            return sequence.reverse_complement()
        return sequence

    def run(self, context_seqs_file, cis_near_distance, genome_fasta, context_seq_len):
        """ Run the annotation pipeline """
        # for performance reasons (mainly seq grepping), results cannot be written on the fly and everything will therefore be stored in a list of lists
        header = [
                "FGID", "Fusion_Gene", "Breakpoint1", "Breakpoint2", "FTID", "context_sequence_id", "context_sequence_100_id", "type",
                "exon_nr", "exon_starts", "exon_ends", "exon_boundary1", "exon_boundary2", "exon_boundary", "bp1_frame", "bp2_frame",
                "frame", "context_sequence", "context_sequence_bp", "neo_peptide_sequence", "neo_peptide_sequence_bp",

#                "context_sequence_wt1", "context_sequence_wt2", "context_sequence_100", "bp1_chr", "bp1_pos", "bp1_strand", "bp2_chr", "bp2_pos", "bp2_strand",
#                "cds_boundary1", "cds_boundary2", "exon_list1", "exon_list2", "cds_list1", "cds_list2", "ft_cds_list1", "ft_cds_list2", "ft_exon_list1",
#                "ft_exon_list2", "seq_list1", "seq_list2", "ft_seq_list1", "ft_seq_list2", "ft_exon_seq_list1", "ft_exon_seq_list2", "wt1_trans_seq",
#                "wt2_trans_seq", "bp1_fusion_trans_seq", "bp2_fusion_trans_seq", "bp1_fusion_trans_seq_ex", "bp2_fusion_trans_seq_ex", "wt1_aa_seq",
#                "wt2_aa_seq", "fusion_nt_seq", "fusion_aa_seq", "is_good_trans1", "is_good_trans2", "trans_biotype1", "trans_biotype2", "gene_biotype1",
#                "gene_biotype2", "description1", "description2", "bp1_frame_at_start", "bp2_frame_at_start", "TSL1", "TSL2", "wt1_cds_no", "wt1_exon_no",
#                "bp1_cds_no", "bp1_exon_no", "wt2_cds_no", "wt2_exon_no", "bp2_cds_no", "bp2_exon_no"
#                ]
#        h2 = [
                "context_sequence_wt1", "context_sequence_wt2", "context_sequence_wt1_bp", "context_sequence_wt2_bp", "context_sequence_100", "bp1_chr",
                "bp1_pos", "bp1_strand", "bp2_chr", "bp2_pos", "bp2_strand", "cds_boundary1", "cds_boundary2", "wt1_exon_pos", "wt2_exon_pos", "ft1_exon_pos",
                "ft2_exon_pos", "wt1_cds_pos", "wt2_cds_pos", "ft1_cds_pos", "ft2_cds_pos", "wt1_exon_seqs", "wt2_exon_seqs", "ft1_exon_seqs", "ft2_exon_seqs",
                "wt1_cds_seqs", "wt2_cds_seqs", "ft1_cds_seqs", "ft2_cds_seqs", "wt1_exon_transcripts", "wt2_exon_transcripts", "ft1_exon_transcripts",
                "ft2_exon_transcripts", "wt1_cds_transcripts", "wt2_cds_transcripts", "ft1_cds_transcripts", "ft2_cds_transcripts", "wt1_peptide",
                "wt2_peptide", "fusion_transcript", "fusion_peptide", "wt1_is_good_transcript", "wt2_is_good_transcript", "wt1_trans_biotype",
                "wt2_trans_biotype", "wt1_gene_biotype", "wt2_gene_biotype", "wt1_description", "wt2_description", "wt1_frame_at_start", "wt2_frame_at_start",
                "wt1_TSL", "wt2_TSL", "wt1_exon_no", "wt2_exon_no", "ft1_exon_no", "ft2_exon_no", "wt1_cds_no", "wt2_cds_no", "ft1_cds_no", "ft2_cds_no"
                ]

        header_idx = range(len(header))
        header_dict = dict(zip(header, header_idx))
        # information in the header is complete and useful for debugging, information in the short_header is identical to the format of the previous
        # version of the fusion annotation and therefore used for downstream processing
        short_header = header[0:21]
        results_lists = []

        # get features overlapping with the bp from the db
        for fgid in self.bp_dict:
            print("Annotating \"{}\"...".format(fgid))
            # create fusion_gene from fgid
            fusion_gene = self.fgid2fusiongene(fgid)
            # split breakpoint info
            bp1 = self.bp_dict[fgid][0]
            bp2 = self.bp_dict[fgid][1]
            bp1_chr, bp1_pos, bp1_strand = bp1.split(":")
            bp2_chr, bp2_pos, bp2_strand = bp2.split(":")
            bp1_pos = int(bp1_pos)
            bp2_pos = int(bp2_pos)
            # define fusion type based on relative breakpoint positions
            fusion_type = self.define_type(cis_near_distance, bp1_chr, bp1_pos, bp1_strand, bp2_chr, bp2_pos, bp2_strand)
            # query db to get exons and matched cds (if available) overlapping the bp as gffutils feature objects
            bp1_exons, bp1_cds = self.get_bp_overlapping_features(bp1_chr, bp1_pos)
            bp2_exons, bp2_cds = self.get_bp_overlapping_features(bp2_chr, bp2_pos)
            # consider all possible combination between transcript1 to transcript2 fusions
            for bp1_efeature, bp1_cfeature in zip(bp1_exons, bp1_cds):
                # check whether breakpoint1 is located on an exon boundary
                exon_boundary1 = self.get_boundary(bp1_pos, bp1_efeature)
                # check whether breakpoint1 is located on a CDS boundary
                cds_boundary1 = self.get_boundary(bp1_pos, bp1_cfeature)
                # get id and additional infos from exon parents transcript and gene
                wt1_trans_id, wt1_trans_biotype, wt1_gene_id, wt1_gene_name, wt1_gene_biotype, wt1_description = self.get_parents(bp1_efeature)
                # create a set to collect "suspicious findings" for transcripts
                if not wt1_trans_id in self.trans_flags_dict:
                    self.trans_flags_dict[wt1_trans_id] = set()
                # get the complete set of exons, corresponding cds (cds at the same index are NA if there is none in the exon) and start frames of the cds
                wt1_exon_pos_list, wt1_cds_pos_list, wt1_cds_frame_list = self.get_wt_codings(wt1_trans_id)
                # get the frame at the start of the first cds and at the breakpoint
                wt1_frame_at_start, wt1_frame_at_bp = self.get_frame(bp1_pos, wt1_cds_pos_list, wt1_cds_frame_list, bp1_strand)
                if wt1_frame_at_start > 0:
                    # flag the transcript if its frame at the start of the first cds is > 0
                    self.trans_flags_dict[wt1_trans_id].add("cds start > 0")
                # remove "NA's" from the cds list
                wt1_cds_pos_list = [x for x in wt1_cds_pos_list if not x == "NA"]
                if len(wt1_cds_pos_list) == 1:
                    # flag the transcript if it contains a single cds
                    self.trans_flags_dict[wt1_trans_id].add("Only 1 cds in bp1")
                # add pre-loaded tsl info
                wt1_tsl = self.get_tsl(wt1_trans_id.split(":")[1])

                # do the same stuff for transcript 2 and (if applicable) combine the information from both transcripts
                for bp2_efeature, bp2_cfeature in zip(bp2_exons, bp2_cds):
                    exon_boundary2 = self.get_boundary(bp2_pos, bp2_efeature)
                    # combined exon boundary information
                    exon_boundary = self.get_boundary2(exon_boundary1, exon_boundary2, bp1_strand, bp2_strand)
                    cds_boundary2 = self.get_boundary(bp2_pos, bp2_cfeature)
                    wt2_trans_id, wt2_trans_biotype, wt2_gene_id, wt2_gene_name, wt2_gene_biotype, wt2_description = self.get_parents(bp2_efeature)
                    if not wt2_trans_id in self.trans_flags_dict:
                        self.trans_flags_dict[wt2_trans_id] = set()
                    # generate ftid from gene names, breakpoints and transcript ids
                    ftid = "_".join([wt1_gene_name, bp1, wt1_trans_id.split(":")[1], wt2_gene_name, bp2, wt2_trans_id.split(":")[1]])
                    wt2_exon_pos_list, wt2_cds_pos_list, wt2_cds_frame_list = self.get_wt_codings(wt2_trans_id)
                    # check whether loci of wt1/wt2 overlap
                    if self.check_exon_overlap(wt1_exon_pos_list, wt2_exon_pos_list, fusion_type):
                        self.trans_flags_dict[wt1_trans_id].add("wt1/wt2 exon overlap")
                        self.trans_flags_dict[wt2_trans_id].add("wt1/wt2 exon overlap")
                    wt2_frame_at_start, wt2_frame_at_bp = self.get_frame(bp2_pos, wt2_cds_pos_list, wt2_cds_frame_list, bp2_strand)
                    if wt2_frame_at_start > 0:
                        self.trans_flags_dict[wt2_trans_id].add("cds start > 0")
                    wt2_cds_pos_list = [x for x in wt2_cds_pos_list if not x == "NA"]
                    wt2_tsl = self.get_tsl(wt2_trans_id.split(":")[1])
                    # combined frame at breakpoint information
                    frame = self.get_frame2(wt1_frame_at_bp, wt2_frame_at_bp)

                    # Get lists of position tuples (cds and exon) which are relevant for the fusion
                    ft1_cds_pos_list, ft2_cds_pos_list = self.get_fusion_feature_list(bp1_pos, bp2_pos, bp1_strand, bp2_strand, wt1_cds_pos_list, wt2_cds_pos_list)
                    ft1_exon_pos_list, ft2_exon_pos_list = self.get_fusion_feature_list(bp1_pos, bp2_pos, bp1_strand, bp2_strand, wt1_exon_pos_list, wt2_exon_pos_list)
                    # Collect positions for which we need to grep sequences
                    self.fill_seq_lookup_dict(bp1_chr, wt1_exon_pos_list)
                    self.fill_seq_lookup_dict(bp2_chr, wt2_exon_pos_list)
                    self.fill_seq_lookup_dict(bp1_chr, ft1_exon_pos_list)
                    self.fill_seq_lookup_dict(bp2_chr, ft2_exon_pos_list)
                    self.fill_seq_lookup_dict(bp1_chr, wt1_cds_pos_list)
                    self.fill_seq_lookup_dict(bp2_chr, wt2_cds_pos_list)
                    self.fill_seq_lookup_dict(bp1_chr, ft1_cds_pos_list)
                    self.fill_seq_lookup_dict(bp2_chr, ft2_cds_pos_list)

                    if frame == "in_frame":
                        exon_nr = len(ft1_cds_pos_list) + len(ft2_cds_pos_list)
                    else:
                        exon_nr = len(ft1_cds_pos_list) + len(ft2_exon_pos_list)
#                    exon_starts = "{}*{}".format("|".join([str(x) for x,y in bp1_cds_pos_list]), "|".join([str(x) for x,y in bp2_cds_pos_list])) # actually cds starts
#                    exon_ends = "{}*{}".format("|".join([str(y) for x,y in bp1_cds_pos_list]), "|".join([str(y) for x,y in bp2_cds_pos_list])) # actually cds ends
                    wt1_is_good_transcript = self.trans_flags_dict[wt1_trans_id]
                    wt2_is_good_transcript = self.trans_flags_dict[wt2_trans_id]
                    results_lists.append([
                            fgid, fusion_gene, bp1, bp2, ftid, "context_sequence_id", "context_sequence_100_id", fusion_type,
                            str(exon_nr), "exon_starts", "exon_ends", exon_boundary1, exon_boundary2, exon_boundary, str(wt1_frame_at_bp), str(wt2_frame_at_bp),
                            frame, "context_sequence", "context_sequence_bp", "neo_peptide_sequence", "neo_peptide_sequence_bp",

                            "context_sequence_wt1", "context_sequence_wt2", "context_sequence_wt1_bp", "context_sequence_wt2_bp", "context_sequence_100", bp1_chr,
                            bp1_pos, bp1_strand, bp2_chr, bp2_pos, bp2_strand, cds_boundary1, cds_boundary2, wt1_exon_pos_list, wt2_exon_pos_list, ft1_exon_pos_list,
                            ft2_exon_pos_list, wt1_cds_pos_list, wt2_cds_pos_list, ft1_cds_pos_list, ft2_cds_pos_list, [], [], [], [],
                            [], [], [], [], "wt1_exon_transcripts", "wt2_exon_transcripts", "ft1_exon_transcripts",
                            "ft2_exon_transcripts", "wt1_cds_transcripts", "wt2_cds_transcripts", "ft1_cds_transcripts", "ft2_cds_transcripts", "wt1_peptide",
                            "wt2_peptide", "fusion_transcript", "fusion_peptide", wt1_is_good_transcript, wt2_is_good_transcript, wt1_trans_biotype,
                            wt2_trans_biotype, wt1_gene_biotype, wt2_gene_biotype, wt1_description, wt2_description, wt1_frame_at_start, wt2_frame_at_start,
                            wt1_tsl, wt2_tsl, len(wt1_exon_pos_list), len(wt2_exon_pos_list), len(ft1_exon_pos_list), len(ft2_exon_pos_list), len(wt1_cds_pos_list), len(wt2_cds_pos_list), len(ft1_cds_pos_list), len(ft2_cds_pos_list)
                            ])

        # grep chromosomes and extract cds sequences from the @cds_seq_dict
        # urla note: the max(..) check in the slicing is probably superfluous as CDS annotations should not start at chromosome position 0
        #            (which would become -1 in @get_frame()). Nevertheless, negative values would completely mess up the slicing...
        #            An additinal check at the end is not required as slicing stops anyway at the end of sequence
        for rec in SeqIO.parse(genome_fasta, "fasta"):
            if rec.id in self.cds_seq_dict:
                print("Grepping candidate regions from chromosome {}".format(rec.id))
                self.cds_seq_dict[rec.id][1] = map(lambda x: rec.seq[max(0, x[0]-1):x[1]], self.cds_seq_dict[rec.id][0])

        # based on the gathered information, create and add sequences to the table
        # urla note: we could have used a pandas df at create functions to apply on the columns to do the following,
        #            but as the following is mainly string manipulation and simple lookups it should not make a big difference
        #            in terms of speed and readability
        last_chr = "" # we save the last chromosome string to avoid some superfluous re-calculations
        count_context_seqs = 0 # the number of context seq records for correct setting of star genome generate parameter
        summed_context_seqs_len = 0 # the cummulative length of context seq records for correct setting of star genome generate parameter
        # we need 3 file writer, 1st for the context seqs, 2nd for the full length transcripts and 3rd for the full length peptides (the last two only for debugging)
        context_fa = "{}.fasta".format(context_seqs_file)
        transcript_fa = "{}_transcript.fasta".format(context_seqs_file)
        peptide_fa = "{}_peptide.fasta".format(context_seqs_file)
        with open(context_fa, "w") as cfasta, open(transcript_fa, "w") as tfasta, open(peptide_fa, "w") as pfasta:
            for result_list in results_lists:
                # create a tempory dict with the cds_pos tuple as keys and the respective genomic sequence as value
                # the dict is restricted to the chromosome of the fusion partner and only updated if the chromosome changes
                if not last_chr == result_list[header_dict["bp1_chr"]]:
                    tmp_dict = dict(zip(self.cds_seq_dict[result_list[header_dict["bp1_chr"]]][0], self.cds_seq_dict[result_list[header_dict["bp1_chr"]]][1]))
                    last_chr = result_list[header_dict["bp1_chr"]]
                # add sequences as gffutils Seq objects to the table
                result_list[header_dict["wt1_exon_seqs"]] = [tmp_dict[feature] for feature in result_list[header_dict["wt1_exon_pos"]]]
                result_list[header_dict["ft1_exon_seqs"]] = [tmp_dict[feature] for feature in result_list[header_dict["ft1_exon_pos"]]]
                result_list[header_dict["wt1_cds_seqs"]] = [tmp_dict[feature] for feature in result_list[header_dict["wt1_cds_pos"]]]
                result_list[header_dict["ft1_cds_seqs"]] = [tmp_dict[feature] for feature in result_list[header_dict["ft1_cds_pos"]]]
                # and concatenate them
                result_list[header_dict["wt1_exon_transcripts"]] = sum(result_list[header_dict["wt1_exon_seqs"]], Seq("", generic_dna))
                result_list[header_dict["ft1_exon_transcripts"]] = sum(result_list[header_dict["ft1_exon_seqs"]], Seq("", generic_dna))
                result_list[header_dict["wt1_cds_transcripts"]] = sum(result_list[header_dict["wt1_cds_seqs"]], Seq("", generic_dna))
                result_list[header_dict["ft1_cds_transcripts"]] = sum(result_list[header_dict["ft1_cds_seqs"]], Seq("", generic_dna))
                # create the first part of the context and fusion sequence and translate wt1 seqs
                translation_shift1 = result_list[header_dict["wt1_frame_at_start"]] # for easier visualization, store in separate var
                trans_table = 1
                if last_chr == "MT":
                    trans_table = 2
                result_list[header_dict["context_sequence"]] = self.get_stranded_seq(result_list[header_dict["ft1_exon_transcripts"]], result_list[header_dict["bp1_strand"]])
                result_list[header_dict["context_sequence_wt1"]] = self.get_stranded_seq(result_list[header_dict["wt1_exon_transcripts"]], result_list[header_dict["bp1_strand"]])
                result_list[header_dict["fusion_transcript"]] = self.get_stranded_seq(result_list[header_dict["ft1_cds_transcripts"]], result_list[header_dict["bp1_strand"]])
                result_list[header_dict["wt1_peptide"]] = self.get_stranded_seq(result_list[header_dict["wt1_cds_transcripts"]], result_list[header_dict["bp1_strand"]])[translation_shift1:].translate(table = trans_table)

                # do the same for fusion partner 2
                if not last_chr == result_list[header_dict["bp2_chr"]]:
                    tmp_dict = dict(zip(self.cds_seq_dict[result_list[header_dict["bp2_chr"]]][0], self.cds_seq_dict[result_list[header_dict["bp2_chr"]]][1]))
                    last_chr = result_list[header_dict["bp2_chr"]]
                result_list[header_dict["wt2_exon_seqs"]] = [tmp_dict[feature] for feature in result_list[header_dict["wt2_exon_pos"]]]
                result_list[header_dict["ft2_exon_seqs"]] = [tmp_dict[feature] for feature in result_list[header_dict["ft2_exon_pos"]]]
                result_list[header_dict["wt2_cds_seqs"]] = [tmp_dict[feature] for feature in result_list[header_dict["wt2_cds_pos"]]]
                result_list[header_dict["ft2_cds_seqs"]] = [tmp_dict[feature] for feature in result_list[header_dict["ft2_cds_pos"]]]
                result_list[header_dict["wt2_exon_transcripts"]] = sum(result_list[header_dict["wt2_exon_seqs"]], Seq("", generic_dna))
                result_list[header_dict["ft2_exon_transcripts"]] = sum(result_list[header_dict["ft2_exon_seqs"]], Seq("", generic_dna))
                result_list[header_dict["wt2_cds_transcripts"]] = sum(result_list[header_dict["wt2_cds_seqs"]], Seq("", generic_dna))
                result_list[header_dict["ft2_cds_transcripts"]] = sum(result_list[header_dict["ft2_cds_seqs"]], Seq("", generic_dna))
                translation_shift2 = result_list[header_dict["wt2_frame_at_start"]]
                # append the second part of the context and fusion sequence and translate wt2 seqs
                trans_table = 1
                if last_chr == "MT":
                    trans_table = 2
                result_list[header_dict["context_sequence"]] += self.get_stranded_seq(result_list[header_dict["ft2_exon_transcripts"]], result_list[header_dict["bp2_strand"]])
                result_list[header_dict["context_sequence_wt2"]] = self.get_stranded_seq(result_list[header_dict["wt2_exon_transcripts"]], result_list[header_dict["bp2_strand"]])
                # for in frame fusions, use the cds sequences, for everything else, use the exons
                if result_list[header_dict["frame"]] == "in_frame":
                    result_list[header_dict["fusion_transcript"]] += self.get_stranded_seq(result_list[header_dict["ft2_cds_transcripts"]], result_list[header_dict["bp2_strand"]])
                else:
                    result_list[header_dict["fusion_transcript"]] += self.get_stranded_seq(result_list[header_dict["ft2_exon_transcripts"]], result_list[header_dict["bp2_strand"]])
                result_list[header_dict["wt2_peptide"]] = self.get_stranded_seq(result_list[header_dict["wt2_cds_transcripts"]], result_list[header_dict["bp2_strand"]])[translation_shift2:].translate(table = trans_table)

                # bp pos in the fusion
                bp_in_fusion_nt = len(result_list[header_dict["ft1_cds_transcripts"]]) # the fusion starts in the fusion transcript with the end of the first part
                bp_in_fusion_nt_ex = len(result_list[header_dict["ft1_exon_transcripts"]]) # the fusion starts in the fusion transcript with the end of the first part
                bp_in_wt1_nt_ex = bp_in_fusion_nt_ex # the breakpoint in wt1 must be identical to the breakpoint in the fusion
                bp_in_wt2_nt_ex = len(result_list[header_dict["wt2_exon_transcripts"]]) - len(result_list[header_dict["ft2_exon_transcripts"]]) # the breakpoint in wt2 is at the position where the ft2 transcripts start
                bp_in_fusion_aa = ((bp_in_fusion_nt - translation_shift1)* 10 // 3) / 10 # the respective bp pos in the aa seq, whereas x.0 / x.3 / x.6 indicate that the breakpoint is at the beginning / after first base / after second base of the underlying codon

                # translate and trim the fusion sequence around the breakpoint
                result_list[header_dict["fusion_peptide"]] = result_list[header_dict["fusion_transcript"]][translation_shift1:].translate(table = 1, to_stop=True)
                result_list[header_dict["context_sequence_100"]] = result_list[header_dict["context_sequence"]][max(0, bp_in_fusion_nt_ex - 100):(bp_in_fusion_nt_ex + 100)]
                result_list[header_dict["context_sequence"]] = result_list[header_dict["context_sequence"]][max(0, bp_in_fusion_nt_ex - context_seq_len):(bp_in_fusion_nt_ex + context_seq_len)]
                result_list[header_dict["context_sequence_wt1"]] = result_list[header_dict["context_sequence_wt1"]][max(0, bp_in_wt1_nt_ex - context_seq_len):(bp_in_wt1_nt_ex + context_seq_len)]
                result_list[header_dict["context_sequence_wt2"]] = result_list[header_dict["context_sequence_wt2"]][max(0, bp_in_wt2_nt_ex - context_seq_len):(bp_in_wt2_nt_ex + context_seq_len)]
                result_list[header_dict["context_sequence_bp"]] = min(context_seq_len, bp_in_fusion_nt_ex)
                result_list[header_dict["context_sequence_wt1_bp"]] = min(context_seq_len, bp_in_wt1_nt_ex)
                result_list[header_dict["context_sequence_wt2_bp"]] = min(context_seq_len, bp_in_wt2_nt_ex)

                result_list[header_dict["neo_peptide_sequence"]] = result_list[header_dict["fusion_peptide"]][max(0, int(bp_in_fusion_aa) - 13):]
                if (int(bp_in_fusion_aa) - 13) < 0:
                    result_list[header_dict["neo_peptide_sequence_bp"]] = bp_in_fusion_aa
                else:
                    result_list[header_dict["neo_peptide_sequence_bp"]] = bp_in_fusion_aa - (int(bp_in_fusion_aa) - 13)
                if result_list[header_dict["frame"]] == "in_frame":
                    result_list[header_dict["neo_peptide_sequence"]] = result_list[header_dict["neo_peptide_sequence"]][:(int(result_list[header_dict["neo_peptide_sequence_bp"]]) + 13)]

                result_list[header_dict["context_sequence_id"]] = xxhash.xxh64(str(result_list[header_dict["context_sequence"]])).hexdigest()
                result_list[header_dict["context_sequence_100_id"]] = xxhash.xxh64(str(result_list[header_dict["context_sequence_100"]])).hexdigest()

                # set exon starts/ends
                # urla: needs revision and simplification, but is working...
                neo_pep_nuc_length = len(result_list[header_dict["neo_peptide_sequence"]]) * 3
                neo_pep_until_bp_nuc = int(round(result_list[header_dict["neo_peptide_sequence_bp"]] * 3, 0))
                if result_list[header_dict["bp1_strand"]] == "+":
                    summed_len = 0
                    for exon_pos, exon_length in enumerate([y-x for x, y in result_list[header_dict["ft1_exon_pos"]]][::-1], 1):
                        summed_len += exon_length
                        if neo_pep_until_bp_nuc <= summed_len:
                            break
                    lookup_exons1 = result_list[header_dict["ft1_exon_pos"]][-exon_pos:]
                else:
                    summed_len = 0
                    for exon_pos, exon_length in enumerate([y-x for x, y in result_list[header_dict["ft1_exon_pos"]]], 1):
                        summed_len += exon_length
                        if neo_pep_until_bp_nuc <= summed_len:
                            break
                    lookup_exons1 = result_list[header_dict["ft1_exon_pos"]][:exon_pos]

                if result_list[header_dict["bp2_strand"]] == "+":
                    summed_len = 0
                    for exon_pos, exon_length in enumerate([y-x for x, y in result_list[header_dict["ft2_exon_pos"]]], 1):
                        summed_len += exon_length
                        if (neo_pep_nuc_length - neo_pep_until_bp_nuc) <= summed_len:
                            break
                    lookup_exons2 = result_list[header_dict["ft2_exon_pos"]][:exon_pos]
                else:
                    summed_len = 0
                    for exon_pos, exon_length in enumerate([y-x for x, y in result_list[header_dict["ft2_exon_pos"]]][::-1], 1):
                        summed_len += exon_length
                        if (neo_pep_nuc_length - neo_pep_until_bp_nuc) <= summed_len:
                            break
                    lookup_exons2 = result_list[header_dict["ft2_exon_pos"]][-exon_pos:]

                result_list[header_dict["exon_starts"]] = "{}*{}".format("|".join([str(x) for x,y in lookup_exons1]), "|".join([str(x) for x,y in lookup_exons2]))
                result_list[header_dict["exon_ends"]] = "{}*{}".format("|".join([str(y) for x,y in lookup_exons1]), "|".join([str(y) for x,y in lookup_exons2]))

                # count context data for star
                # we perform tsl exclusion for the context seqs, but keep everything for the debugging data
                if not (result_list[header_dict["wt1_TSL"]] in ["4", "5"] or result_list[header_dict["wt2_TSL"]] in ["4", "5"]):
                    count_context_seqs += 3
                    summed_context_seqs_len += sum(map(len, [
                            result_list[header_dict["context_sequence"]],
                            result_list[header_dict["context_sequence_wt1"]],
                            result_list[header_dict["context_sequence_wt2"]]
                            ]))
                    # write sequences for requantification
                    SeqIO.write(SeqRecord(result_list[header_dict["context_sequence"]], id="{}_{}_{}_ft".format(result_list[header_dict["FTID"]], result_list[header_dict["context_sequence_id"]], result_list[header_dict["context_sequence_bp"]]), name="", description=""), cfasta, "fasta")
                    SeqIO.write(SeqRecord(result_list[header_dict["context_sequence_wt1"]], id="{}_{}_{}_wt1".format(result_list[header_dict["FTID"]], result_list[header_dict["context_sequence_id"]], result_list[header_dict["context_sequence_wt1_bp"]]), name="", description=""), cfasta, "fasta")
                    SeqIO.write(SeqRecord(result_list[header_dict["context_sequence_wt2"]], id="{}_{}_{}_wt2".format(result_list[header_dict["FTID"]], result_list[header_dict["context_sequence_id"]], result_list[header_dict["context_sequence_wt2_bp"]]), name="", description=""), cfasta, "fasta")

                SeqIO.write(SeqRecord(result_list[header_dict["fusion_transcript"]], id="{}_ft_{}".format(result_list[header_dict["FTID"]], bp_in_fusion_nt), name="", description=""), tfasta, "fasta")
                SeqIO.write(SeqRecord(result_list[header_dict["fusion_peptide"]], id="{}_{}_{}".format(result_list[header_dict["FTID"]], bp_in_fusion_aa, result_list[header_dict["type"]]), name="", description=""), pfasta, "fasta")

            # write sequences to file

            # experimental
            if not len(result_list[header_dict["wt1_cds_transcripts"]]) % 3 == 0:
                result_list[header_dict["wt1_is_good_transcript"]].add("wt1 seq % 3 != 0")
            if not len(result_list[header_dict["wt2_cds_transcripts"]]) % 3 == 0:
                result_list[header_dict["wt2_is_good_transcript"]].add("wt2 seq % 3 != 0")

            result_list[header_dict["wt1_cds_transcripts"]] = str(result_list[header_dict["wt1_cds_transcripts"]])
            result_list[header_dict["wt2_cds_transcripts"]] = str(result_list[header_dict["wt2_cds_transcripts"]])

        # write everything to file
        # for easy debugging, we will write two files. The first contains all informations, the second is filtered and contains only downstream required columns
        with open("{}.debug".format(context_seqs_file), "w") as outfile:
            outfile.write("{}\n".format(";".join(header)))
            for result_line in results_lists:
                outfile.write("{}\n".format(";".join(map(str, result_line))))
        with open(context_seqs_file, "w") as outfile:
            outfile.write("{}\n".format(";".join(short_header)))
            for result_line in results_lists:
                # filter all fusions with transcripts which have a TSL of 4 or 5 (keeping 1-3, 6(liftover), NA)
                if result_line[header_dict["wt1_TSL"]] in ["4", "5"] or result_line[header_dict["wt2_TSL"]] in ["4", "5"]:
                    continue
                outfile.write("{}\n".format(";".join([str(result_line[header_dict[col_name]]) for col_name in short_header])))
        # write fasta.info file for star genome generate
        with open("{}.fasta.info".format(context_seqs_file), "w") as ifasta:
            ifasta.write("{}\n".format(count_context_seqs))
            ifasta.write("{}\n".format(summed_context_seqs_len))


def main():
    """Parse command line arguments and start script"""
    parser = argparse.ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('--detected_fusions', dest='detected_fusions', help='detected_fusions.csv file from easyfuse', required=True)
    parser.add_argument('--annotation_db', dest='annotation_db', help='gff3 derived annotation database from gffutils', required=True)
    parser.add_argument('--out_csv', dest='out_csv', help='context_seqs.csv output for easyfuse', required=True)
    parser.add_argument('--genome_fasta', dest='genome_fasta', help='Genome fasta file for sequence grepping', required=True)
    parser.add_argument('--tsl_info', dest='tsl_info', help='File containing transcript support level (TSL). Should have been generated by gtf2tsl.py', required=True)
    parser.add_argument('--cis_near_dist', dest='cis_near_dist', help='Genomic distance after which a cis_near becomes a cis_far fusion', required=True)
    parser.add_argument('--context_seq_length', dest='context_seq_len', help='The maximal length of sequence context (exon) at either site of the bp', required=True)
    args = parser.parse_args()

#    bp_in_file = "/mnt/bfx/urla/easyfuse_validationData/easyfuse_standalone_urla/Sample_EasyValiS01TR_NFC15_Rep1/fetchdata/fd_1_tool/fetched_fusions/Detected_Fusions.csv"
#    db_in = "/mnt/bfx/IVAC_DEV/src/FUS/GenomeData/Ensembl_GRCh37_v90/Homo_sapiens.GRCh37.87.annot.db"
#    out_file = "/mnt/bfx/urla/easyfuse_validationData/easyfuse_standalone_urla/newAnnotTest/Context_Seqs_urla_S01_R1.csv"
#    genome_fasta = "/mnt/bfx/IVAC_DEV/src/FUS/GenomeData/Ensembl_GRCh37_v90/Homo_sapiens.GRCh37.dna.primary_assembly.fasta"
#    tsl_in = "/mnt/bfx/IVAC_DEV/src/FUS/GenomeData/Ensembl_GRCh37_v90/Homo_sapiens.GRCh37.87.tsl_info"

    fusannot = FusionAnnotation(args.annotation_db, args.detected_fusions, args.tsl_info)
    fusannot.run(args.out_csv, int(args.cis_near_dist), args.genome_fasta, int(args.context_seq_len))

if __name__ == '__main__':
    main()
