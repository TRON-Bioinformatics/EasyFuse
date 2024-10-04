#!/usr/bin/env python3

"""
Annotate predicted fusion genes. This involves:
    1) Identify transcripts possibly involved in the fusion
    2) Determine the translation frame at the breakpoint position
    3) Reconstruct possible fusion transcripts as well as the respective wildtype background
    4) Translate transcripts
    5) Select regions of interest for neopeptide searches
    6) Collect information on the kind of fusion which might be relevant for downstream
       selection processes
    7) Filter likely artifacts based on the quality of the annotated transcripts
@author: TRON (PASO), BNT (URLA)
@version: 241004
"""

from argparse import ArgumentParser

# pylint: disable=E0401
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
from logzero import logger
import xxhash


from file_headers import FULL_ANNOTATION_HEADER
from io_methods import breakpoints_to_dict
from io_methods import load_tsl_data
from boundary_validation import get_boundary, get_boundary2
from exon_validation import check_exon_overlap
from exon_validation import get_wt_codings
from exon_validation import get_bp_overlapping_features
from exon_validation import get_parents
from frame_validation import get_frame, get_frame2
from fusion_validation import define_type, get_fusion_feature_list
from sequence_validation import get_stranded_seq



class FusionAnnotation:
    """Annotation of predicted fusion genes soley based on the breakpoint information"""

    def __init__(self, db_file, bp_file, tsl_file):
        """Parameter initialization"""
        self.bp_dict = breakpoints_to_dict(bp_file)
        self.db = gffutils.FeatureDB(db_file)
        self.tsl_dict = load_tsl_data(tsl_file)
        self.cds_seq_dict = {}
        self.trans_flags_dict = {}


    def get_tsl(self, trans_id):
        """Get the transcript support level (tsl) from the pre-loaded dict"""
        # the tsl is not available for grch37 data (at least not v87/90)
        # and I'm therefore using the tsl info from grch38 v90
        # With the newer data, however, it is possible that some transcript ids
        # are not in the database (because the respective
        # transcript annotations were deleted or modified significantly).
        # Although the "correct" tsl for a missing trans_id would be "NA",
        # I'm using "6" in order to allow differentiation from
        # existing, "normal" tsl's which range from 1-5 and are NA
        # for transcripts of unknown support
        if not trans_id in self.tsl_dict:
            return "6"
        return self.tsl_dict[trans_id]


    def fill_seq_lookup_dict(self, chrom, cds_list):
        """
        Create/fill a dictionary of genomic position tuple with chromosomes as keys 
        and a list of two lists as values. Genomic coordinates are loaded into the first list 
        and the second one is later complemented by the respective sequence
        """
        if not chrom in self.cds_seq_dict:
            self.cds_seq_dict[chrom] = [[], []]
        for cds_coord in cds_list:
            if not cds_coord in self.cds_seq_dict[chrom][0]:
                self.cds_seq_dict[chrom][0].append(cds_coord)


    def run(
        self,
        context_seqs_file,
        cis_near_distance,
        genome_fasta,
        context_seq_len,
        tsl_filter_level,
    ):
        """Run the annotation pipeline"""
        # for performance reasons (mainly seq grepping), results cannot be written
        # on the fly and everything will therefore be stored in a list of lists
        logger.info("Running annotation pipeline")

        header_idx = range(len(FULL_ANNOTATION_HEADER))
        header_dict = dict(zip(FULL_ANNOTATION_HEADER, header_idx))
        # information in the header is complete and useful for debugging,
        # information in the short_header is identical
        # to the format of the previous
        # version of the fusion annotation and therefore used for downstream processing
        short_header = FULL_ANNOTATION_HEADER[0:25]
        results_lists = []

        # get features overlapping with the bp from the db
        for bpid in sorted(self.bp_dict):
            logger.debug("Annotating '%s'...", bpid)
            # split breakpoint info
            bp1 = self.bp_dict[bpid][0]
            bp2 = self.bp_dict[bpid][1]
            bp1_chr, bp1_pos, bp1_strand = bp1.split(":")
            bp2_chr, bp2_pos, bp2_strand = bp2.split(":")
            bp1_pos = int(bp1_pos)
            bp2_pos = int(bp2_pos)
            # define fusion type based on relative breakpoint positions
            fusion_type = define_type(
                cis_near_distance,
                bp1_chr,
                bp1_pos,
                bp1_strand,
                bp2_chr,
                bp2_pos,
                bp2_strand,
            )
            # query db to get exons and matched cds (if available) overlapping
            # the bp as gffutils feature objects
            bp1_exons, bp1_cds = get_bp_overlapping_features(self.db, bp1_chr, bp1_pos)
            bp2_exons, bp2_cds = get_bp_overlapping_features(self.db, bp2_chr, bp2_pos)
            # consider all possible combination between transcript1 to transcript2 fusions
            for bp1_efeature, bp1_cfeature in zip(bp1_exons, bp1_cds):
                # check whether breakpoint1 is located on an exon boundary
                exon_boundary1 = get_boundary(bp1_pos, bp1_efeature)
                # check whether breakpoint1 is located on a CDS boundary
                cds_boundary1 = get_boundary(bp1_pos, bp1_cfeature)
                # get id and additional infos from exon parents transcript and gene
                (
                    wt1_trans_id,
                    wt1_trans_biotype,
                    wt1_gene_name,
                    wt1_gene_biotype,
                    wt1_description,
                ) = get_parents(self.db, bp1_efeature)
                # create a set to collect "suspicious findings" for transcripts
                if wt1_trans_id not in self.trans_flags_dict:
                    self.trans_flags_dict[wt1_trans_id] = set()
                # get the complete set of exons, corresponding cds
                # (cds at the same index are NA if there is none in
                # the exon) and start frames of the cds
                (wt1_exon_pos_list, wt1_cds_pos_list, wt1_cds_frame_list) = get_wt_codings(
                    self.db,
                    self.trans_flags_dict,
                    wt1_trans_id
                )
                # get the frame at the start of the first cds and at the breakpoint
                wt1_frame_at_start, wt1_frame_at_bp = get_frame(
                    bp1_pos, wt1_cds_pos_list, wt1_cds_frame_list, bp1_strand
                )
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

                # do the same stuff for transcript 2 and (if applicable)
                # combine the information from both transcripts
                for bp2_efeature, bp2_cfeature in zip(bp2_exons, bp2_cds):
                    exon_boundary2 = get_boundary(bp2_pos, bp2_efeature)
                    # combined exon boundary information
                    exon_boundary = get_boundary2(
                        exon_boundary1, exon_boundary2, bp1_strand, bp2_strand
                    )
                    cds_boundary2 = get_boundary(bp2_pos, bp2_cfeature)
                    (
                        wt2_trans_id,
                        wt2_trans_biotype,
                        wt2_gene_name,
                        wt2_gene_biotype,
                        wt2_description,
                    ) = get_parents(self.db, bp2_efeature)
                    if wt2_trans_id not in self.trans_flags_dict:
                        self.trans_flags_dict[wt2_trans_id] = set()
                    # generate ftid from gene names, breakpoints and transcript ids
                    ftid = "_".join(
                        [
                            wt1_gene_name,
                            bp1,
                            wt1_trans_id.split(":")[1],
                            wt2_gene_name,
                            bp2,
                            wt2_trans_id.split(":")[1],
                        ]
                    )
                    # change fgid and create fusion_gene based on the easyfuse annotation
                    # fgid_ef = "_".join([wt1_gene_name, bp1, wt2_gene_name, bp2])
                    fusion_gene_ef = "_".join([wt1_gene_name, wt2_gene_name])

                    (wt2_exon_pos_list, wt2_cds_pos_list, wt2_cds_frame_list) = get_wt_codings(
                        self.db,
                        self.trans_flags_dict,
                        wt2_trans_id
                    )
                    # check whether loci of wt1/wt2 overlap
                    if check_exon_overlap(
                        wt1_exon_pos_list, wt2_exon_pos_list, fusion_type
                    ):
                        self.trans_flags_dict[wt1_trans_id].add("wt1/wt2 exon overlap")
                        self.trans_flags_dict[wt2_trans_id].add("wt1/wt2 exon overlap")
                    wt2_frame_at_start, wt2_frame_at_bp = get_frame(
                        bp2_pos, wt2_cds_pos_list, wt2_cds_frame_list, bp2_strand
                    )
                    if wt2_frame_at_start > 0:
                        self.trans_flags_dict[wt2_trans_id].add("cds start > 0")
                    wt2_cds_pos_list = [x for x in wt2_cds_pos_list if not x == "NA"]
                    wt2_tsl = self.get_tsl(wt2_trans_id.split(":")[1])
                    # combined frame at breakpoint information
                    frame = get_frame2(wt1_frame_at_bp, wt2_frame_at_bp)

                    # Get lists of position tuples (cds and exon) which are relevant for the fusion
                    ft1_cds_pos_list, ft2_cds_pos_list = get_fusion_feature_list(
                        bp1_pos,
                        bp2_pos,
                        bp1_strand,
                        bp2_strand,
                        wt1_cds_pos_list,
                        wt2_cds_pos_list,
                    )
                    ft1_exon_pos_list, ft2_exon_pos_list = get_fusion_feature_list(
                        bp1_pos,
                        bp2_pos,
                        bp1_strand,
                        bp2_strand,
                        wt1_exon_pos_list,
                        wt2_exon_pos_list,
                    )
                    # Collect positions for which we need to grep sequences
                    self.fill_seq_lookup_dict(bp1_chr, wt1_exon_pos_list)
                    #self.fill_seq_lookup_dict(bp1_chr, wt1_exon_until_bp_pos_list)
                    self.fill_seq_lookup_dict(bp2_chr, wt2_exon_pos_list)
                    #self.fill_seq_lookup_dict(bp2_chr, wt2_exon_until_bp_pos_list)
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

                    ft1_exon_nr = len(ft1_exon_pos_list)
                    ft2_exon_nr = len(wt1_exon_pos_list) - len(ft2_exon_pos_list)
                    wt1_is_good_transcript = self.trans_flags_dict[wt1_trans_id]
                    wt2_is_good_transcript = self.trans_flags_dict[wt2_trans_id]
                    results_lists.append(
                        [
                            bpid,
                            fusion_gene_ef,
                            bp1,
                            bp2,
                            ftid,
                            "context_sequence_id",
                            "context_sequence_100_id",
                            fusion_type,
                            str(exon_nr),
                            str(ft1_exon_nr),
                            str(ft2_exon_nr),
                            "exon_starts",
                            "exon_ends",
                            exon_boundary1,
                            exon_boundary2,
                            exon_boundary,
                            str(wt1_frame_at_bp),
                            str(wt2_frame_at_bp),
                            frame,
                            "context_sequence",
                            "context_sequence_bp",
                            "neo_peptide_sequence",
                            "neo_peptide_sequence_bp",
                            "fusion_protein_sequence",
                            "fusion_protein_sequence_bp",
                            "context_sequence_wt1",
                            "context_sequence_wt2",
                            "context_sequence_wt1_bp",
                            "context_sequence_wt2_bp",
                            "context_sequence_100",
                            bp1_chr,
                            bp1_pos,
                            bp1_strand,
                            bp2_chr,
                            bp2_pos,
                            bp2_strand,
                            cds_boundary1,
                            cds_boundary2,
                            wt1_exon_pos_list,
                            wt2_exon_pos_list,
                            ft1_exon_pos_list,
                            ft2_exon_pos_list,
                            wt1_cds_pos_list,
                            wt2_cds_pos_list,
                            ft1_cds_pos_list,
                            ft2_cds_pos_list,
                            [],
                            [],
                            [],
                            [],
                            [],
                            [],
                            [],
                            [],
                            "wt1_exon_transcripts",
                            "wt2_exon_transcripts",
                            "ft1_exon_transcripts",
                            "ft2_exon_transcripts",
                            "wt1_cds_transcripts",
                            "wt2_cds_transcripts",
                            "ft1_cds_transcripts",
                            "ft2_cds_transcripts",
                            "wt1_peptide",
                            "wt2_peptide",
                            "fusion_transcript",
                            "fusion_peptide",
                            wt1_is_good_transcript,
                            wt2_is_good_transcript,
                            wt1_trans_biotype,
                            wt2_trans_biotype,
                            wt1_gene_biotype,
                            wt2_gene_biotype,
                            wt1_description,
                            wt2_description,
                            wt1_frame_at_start,
                            wt2_frame_at_start,
                            wt1_tsl,
                            wt2_tsl,
                            len(wt1_exon_pos_list),
                            len(wt2_exon_pos_list),
                            len(ft1_exon_pos_list),
                            len(ft2_exon_pos_list),
                            len(wt1_cds_pos_list),
                            len(wt2_cds_pos_list),
                            len(ft1_cds_pos_list),
                            len(ft2_cds_pos_list),
                            f"{bp1_chr}:{wt1_exon_pos_list[0][0]}:{wt1_exon_pos_list[-1][1]}",
                            f"{bp2_chr}:{wt2_exon_pos_list[0][0]}:{wt2_exon_pos_list[-1][1]}",
                            "annotation_bias",
                        ]
                    )

        # grep chromosomes and extract cds sequences from the @cds_seq_dict
        # urla note: the max(..) check in the slicing is probably superfluous
        # as CDS annotations should not start at chromosome position 0
        # (which would become -1 in @get_frame()). Nevertheless, negative values would
        # completely mess up the slicing...
        # An additinal check at the end is not required as
        # slicing stops anyway at the end of sequence
        for record in SeqIO.parse(genome_fasta, "fasta"):
            if record.id in self.cds_seq_dict:
                logger.info(
                    "Grepping candidate regions from chromosome %s", record.id
                )

                self.cds_seq_dict[record.id][1] = [
                    record.seq[max(0, val[0] - 1) : val[1]]
                    for val in self.cds_seq_dict[record.id][0]
                ]

        # based on the gathered information, create and add sequences to the table
        # urla note: we could have used a pandas df at create functions to
        # apply on the columns to do the following,
        # but as the following is mainly string manipulation and simple lookups
        # should not make a big difference in terms of speed and readability
        last_chr = ""
        # we save the last chromosome string to avoid some superfluous re-calculations
        # due to annotation biases between different prediction tools,
        # it is possible that the same ftid is created twice
        # as this is 100% redundant information, we will skip such records.
        # In order to guarantee that this is working as intended,
        # we introduce the ftid plus which is the ftid plus appended seq hash
        ftid_plus_set = set()
        tmp_dict = None
        for result_list in results_lists:
            # create a tempory dict with the cds_pos tuple as keys and
            # the respective genomic sequence as value
            # the dict is restricted to the chromosome of the fusion partner and
            # only updated if the chromosome changes

            if not last_chr == result_list[header_dict["bp1_chr"]]:
                tmp_dict = dict(
                    zip(
                        self.cds_seq_dict[result_list[header_dict["bp1_chr"]]][0],
                        self.cds_seq_dict[result_list[header_dict["bp1_chr"]]][1],
                    )
                )
                last_chr = result_list[header_dict["bp1_chr"]]
            # add sequences as gffutils Seq objects to the table

            result_list[header_dict["wt1_exon_seqs"]] = []
            for feature in result_list[header_dict["wt1_exon_pos"]]:
                result_list[header_dict["wt1_exon_seqs"]].append(tmp_dict[feature])
            result_list[header_dict["ft1_exon_seqs"]] = [
                tmp_dict[feature]
                for feature in result_list[header_dict["ft1_exon_pos"]]
            ]
            result_list[header_dict["wt1_cds_seqs"]] = [
                tmp_dict[feature] for feature in result_list[header_dict["wt1_cds_pos"]]
            ]
            result_list[header_dict["ft1_cds_seqs"]] = [
                tmp_dict[feature] for feature in result_list[header_dict["ft1_cds_pos"]]
            ]
            # and concatenate them
            result_list[header_dict["wt1_exon_transcripts"]] = sum(
                result_list[header_dict["wt1_exon_seqs"]], Seq("")
            )
            result_list[header_dict["ft1_exon_transcripts"]] = sum(
                result_list[header_dict["ft1_exon_seqs"]], Seq("")
            )
            result_list[header_dict["wt1_cds_transcripts"]] = sum(
                result_list[header_dict["wt1_cds_seqs"]], Seq("")
            )
            result_list[header_dict["ft1_cds_transcripts"]] = sum(
                result_list[header_dict["ft1_cds_seqs"]], Seq("")
            )
            # create the first part of the context and fusion sequence and translate wt1 seqs
            translation_shift1 = result_list[
                header_dict["wt1_frame_at_start"]
            ]  # for easier visualization, store in separate var
            trans_table = 1
            if last_chr == "MT":
                trans_table = 2
            result_list[header_dict["context_sequence"]] = get_stranded_seq(
                result_list[header_dict["ft1_exon_transcripts"]],
                result_list[header_dict["bp1_strand"]],
            )
            result_list[header_dict["context_sequence_wt1"]] = get_stranded_seq(
                result_list[header_dict["wt1_exon_transcripts"]],
                result_list[header_dict["bp1_strand"]],
            )
            result_list[header_dict["fusion_transcript"]] = get_stranded_seq(
                result_list[header_dict["ft1_cds_transcripts"]],
                result_list[header_dict["bp1_strand"]],
            )
            result_list[header_dict["wt1_peptide"]] = get_stranded_seq(
                result_list[header_dict["wt1_cds_transcripts"]],
                result_list[header_dict["bp1_strand"]],
            )[translation_shift1:].translate(table=trans_table)

            # do the same for fusion partner 2
            if not last_chr == result_list[header_dict["bp2_chr"]]:
                tmp_dict = dict(
                    zip(
                        self.cds_seq_dict[result_list[header_dict["bp2_chr"]]][0],
                        self.cds_seq_dict[result_list[header_dict["bp2_chr"]]][1],
                    )
                )
                last_chr = result_list[header_dict["bp2_chr"]]
            result_list[header_dict["wt2_exon_seqs"]] = [
                tmp_dict[feature]
                for feature in result_list[header_dict["wt2_exon_pos"]]
            ]
            result_list[header_dict["ft2_exon_seqs"]] = [
                tmp_dict[feature]
                for feature in result_list[header_dict["ft2_exon_pos"]]
            ]
            result_list[header_dict["wt2_cds_seqs"]] = [
                tmp_dict[feature] for feature in result_list[header_dict["wt2_cds_pos"]]
            ]
            result_list[header_dict["ft2_cds_seqs"]] = [
                tmp_dict[feature] for feature in result_list[header_dict["ft2_cds_pos"]]
            ]
            result_list[header_dict["wt2_exon_transcripts"]] = sum(
                result_list[header_dict["wt2_exon_seqs"]], Seq("")
            )
            result_list[header_dict["ft2_exon_transcripts"]] = sum(
                result_list[header_dict["ft2_exon_seqs"]], Seq("")
            )
            result_list[header_dict["wt2_cds_transcripts"]] = sum(
                result_list[header_dict["wt2_cds_seqs"]], Seq("")
            )
            result_list[header_dict["ft2_cds_transcripts"]] = sum(
                result_list[header_dict["ft2_cds_seqs"]], Seq("")
            )
            translation_shift2 = result_list[header_dict["wt2_frame_at_start"]]
            # append the second part of the context and fusion sequence and translate wt2 seqs
            trans_table = 1
            if last_chr == "MT":
                trans_table = 2
            result_list[header_dict["context_sequence"]] += get_stranded_seq(
                result_list[header_dict["ft2_exon_transcripts"]],
                result_list[header_dict["bp2_strand"]],
            )
            result_list[header_dict["context_sequence_wt2"]] = get_stranded_seq(
                result_list[header_dict["wt2_exon_transcripts"]],
                result_list[header_dict["bp2_strand"]],
            )
            # for in frame fusions, use the cds sequences, for everything else, use the exons
            if result_list[header_dict["frame"]] == "in_frame":
                result_list[header_dict["fusion_transcript"]] += get_stranded_seq(
                    result_list[header_dict["ft2_cds_transcripts"]],
                    result_list[header_dict["bp2_strand"]],
                )
            else:
                result_list[header_dict["fusion_transcript"]] += get_stranded_seq(
                    result_list[header_dict["ft2_exon_transcripts"]],
                    result_list[header_dict["bp2_strand"]],
                )
            result_list[header_dict["wt2_peptide"]] = get_stranded_seq(
                result_list[header_dict["wt2_cds_transcripts"]],
                result_list[header_dict["bp2_strand"]],
            )[translation_shift2:].translate(table=trans_table)

            # bp pos in the fusion
            # the fusion starts in the fusion transcript with the end of the first part
            bp_in_fusion_nt = len(result_list[header_dict["ft1_cds_transcripts"]])
            # the fusion starts in the fusion transcript with the end of the first part
            bp_in_fusion_nt_ex = len(result_list[header_dict["ft1_exon_transcripts"]])
            # the breakpoint in wt1 must be identical to the breakpoint in the fusion
            bp_in_wt1_nt_ex = bp_in_fusion_nt_ex
            # the breakpoint in wt2 is at the position where the ft2 transcripts start
            bp_in_wt2_nt_ex = len(
                result_list[header_dict["wt2_exon_transcripts"]]
            ) - len(result_list[header_dict["ft2_exon_transcripts"]])
            # the respective bp pos in the aa seq, whereas x.0 / x.3 / x.6
            # indicate that the breakpoint is at the beginning / after first base /
            # after second base of the underlying codon
            bp_in_fusion_aa = (
                (bp_in_fusion_nt - translation_shift1) * 10.0 // 3
            ) / 10.0

            # translate and trim the fusion sequence around the breakpoint
            result_list[header_dict["fusion_peptide"]] = result_list[
                header_dict["fusion_transcript"]
            ][translation_shift1:].translate(table=1, to_stop=True)
            result_list[header_dict["context_sequence_100"]] = result_list[
                header_dict["context_sequence"]
            ][max(0, bp_in_fusion_nt_ex - 100) : (bp_in_fusion_nt_ex + 100)]
            result_list[header_dict["context_sequence"]] = result_list[
                header_dict["context_sequence"]
            ][
                max(0, bp_in_fusion_nt_ex - context_seq_len) : (
                    bp_in_fusion_nt_ex + context_seq_len
                )
            ]
            result_list[header_dict["context_sequence_wt1"]] = result_list[
                header_dict["context_sequence_wt1"]
            ][
                max(0, bp_in_wt1_nt_ex - context_seq_len) : (
                    bp_in_wt1_nt_ex + context_seq_len
                )
            ]
            result_list[header_dict["context_sequence_wt2"]] = result_list[
                header_dict["context_sequence_wt2"]
            ][
                max(0, bp_in_wt2_nt_ex - context_seq_len) : (
                    bp_in_wt2_nt_ex + context_seq_len
                )
            ]
            result_list[header_dict["context_sequence_bp"]] = min(
                context_seq_len, bp_in_fusion_nt_ex
            )
            result_list[header_dict["context_sequence_wt1_bp"]] = min(
                context_seq_len, bp_in_wt1_nt_ex
            )
            result_list[header_dict["context_sequence_wt2_bp"]] = min(
                context_seq_len, bp_in_wt2_nt_ex
            )

            result_list[header_dict["neo_peptide_sequence"]] = result_list[
                header_dict["fusion_peptide"]
            ][max(0, int(bp_in_fusion_aa) - 13) :]

            result_list[header_dict["fusion_protein_sequence"]] = result_list[
                header_dict["fusion_peptide"]
            ]

            result_list[header_dict["fusion_protein_sequence_bp"]] = round(
                    bp_in_fusion_aa, 1
                )

            if (int(bp_in_fusion_aa) - 13) < 0:
                result_list[header_dict["neo_peptide_sequence_bp"]] = round(
                    bp_in_fusion_aa, 1
                )
            else:
                result_list[header_dict["neo_peptide_sequence_bp"]] = round(
                    bp_in_fusion_aa - (int(bp_in_fusion_aa) - 13), 1
                )
            if result_list[header_dict["frame"]] == "in_frame":
                result_list[header_dict["neo_peptide_sequence"]] = result_list[
                    header_dict["neo_peptide_sequence"]
                ][: (int(result_list[header_dict["neo_peptide_sequence_bp"]]) + 13)]

            result_list[header_dict["context_sequence_id"]] = xxhash.xxh64(
                str(result_list[header_dict["context_sequence"]])
            ).hexdigest()
            result_list[header_dict["context_sequence_100_id"]] = xxhash.xxh64(
                str(result_list[header_dict["context_sequence_100"]])
            ).hexdigest()

            ftid = result_list[header_dict["FTID"]]
            ctx_id = result_list[header_dict["context_sequence_id"]]
            ftid_plus = f"{ftid}_{ctx_id}"
            if ftid_plus in ftid_plus_set:
                result_list[header_dict["annotation_bias"]] = True
            else:
                result_list[header_dict["annotation_bias"]] = False
                ftid_plus_set.add(ftid_plus)

            # set exon starts/ends
            # urla: needs revision and simplification, but is working...
            neo_pep_nuc_length = (
                len(result_list[header_dict["neo_peptide_sequence"]]) * 3
            )
            neo_pep_until_bp_nuc = int(
                round(result_list[header_dict["neo_peptide_sequence_bp"]] * 3, 0)
            )
            exon_pos = 1
            if result_list[header_dict["bp1_strand"]] == "+":
                summed_len = 0
                for exon_pos, exon_length in enumerate(
                    [y - x for x, y in result_list[header_dict["ft1_exon_pos"]]][::-1],
                    1,
                ):
                    summed_len += exon_length
                    if neo_pep_until_bp_nuc <= summed_len:
                        break
                lookup_exons1 = result_list[header_dict["ft1_exon_pos"]][-exon_pos:]
            else:
                summed_len = 0
                for exon_pos, exon_length in enumerate(
                    [y - x for x, y in result_list[header_dict["ft1_exon_pos"]]], 1
                ):
                    summed_len += exon_length
                    if neo_pep_until_bp_nuc <= summed_len:
                        break
                lookup_exons1 = result_list[header_dict["ft1_exon_pos"]][:exon_pos]

            if result_list[header_dict["bp2_strand"]] == "+":
                summed_len = 0
                for exon_pos, exon_length in enumerate(
                    [y - x for x, y in result_list[header_dict["ft2_exon_pos"]]], 1
                ):
                    summed_len += exon_length
                    if (neo_pep_nuc_length - neo_pep_until_bp_nuc) <= summed_len:
                        break
                lookup_exons2 = result_list[header_dict["ft2_exon_pos"]][:exon_pos]
            else:
                summed_len = 0
                for exon_pos, exon_length in enumerate(
                    [y - x for x, y in result_list[header_dict["ft2_exon_pos"]]][::-1],
                    1,
                ):
                    summed_len += exon_length
                    if (neo_pep_nuc_length - neo_pep_until_bp_nuc) <= summed_len:
                        break
                lookup_exons2 = result_list[header_dict["ft2_exon_pos"]][-exon_pos:]

            exon_starts1_str = "|".join([str(x) for x, y in lookup_exons1])
            exon_starts2_str = "|".join([str(x) for x, y in lookup_exons2])
            result_list[header_dict["exon_starts"]] = f"{exon_starts1_str}*{exon_starts2_str}"

            exon_ends1_str = "|".join([str(y) for x, y in lookup_exons1])
            exon_ends2_str = "|".join([str(y) for x, y in lookup_exons2])
            result_list[header_dict["exon_ends"]] = f"{exon_ends1_str}*{exon_ends2_str}"

            # experimental
            if not len(result_list[header_dict["wt1_cds_transcripts"]]) % 3 == 0:
                result_list[header_dict["wt1_is_good_transcript"]].add(
                    "wt1 seq % 3 != 0"
                )
            if not len(result_list[header_dict["wt2_cds_transcripts"]]) % 3 == 0:
                result_list[header_dict["wt2_is_good_transcript"]].add(
                    "wt2 seq % 3 != 0"
                )

            result_list[header_dict["wt1_cds_transcripts"]] = str(
                result_list[header_dict["wt1_cds_transcripts"]]
            )
            result_list[header_dict["wt2_cds_transcripts"]] = str(
                result_list[header_dict["wt2_cds_transcripts"]]
            )

        # write everything to a couple of files
        # 1) context.csv.debug:
        # all fusions -> used for debugging
        # 2) context.csv:
        # filtered fusions -> contains only downstream required columns
        # 3) context_seqs.csv.fasta:
        # contain the fasta sequences of the ft/wt1/wt2 context sequences,
        # the full length reconstructed transcripts of the gene fusion and
        # the full length translated peptide sequence of the gene fusion
        # 4) context_seqs.csv_transcripts.fasta (for debugging)
        # 5) context_seqs.csv_peptide.fasta (for debugging)
        # 6) context_seqs.csv.fasta.info :
        # Length and count of fasta sequences in context_seqs.csv.fasta

        # the number of context seq records for correct setting of star genome generate parameter
        count_context_seqs = 0
        # the cummulative length of context seq records for correct
        # setting of star genome generate parameter
        summed_context_seqs_len = 0
        with open(f"{context_seqs_file}.debug", "w", encoding="utf8") as outfile_debug, open(
            context_seqs_file, "w", encoding="utf8"
        ) as outfile, open(f"{context_seqs_file}.fasta", "w", encoding="utf8") as cfasta, open(
            f"{context_seqs_file}_transcript.fasta", "w", encoding="utf8"
        ) as tfasta, open(
            f"{context_seqs_file}_peptide.fasta", "w", encoding="utf8"
        ) as pfasta:
            # write header to csv files
            full_annotation_header_str = ";".join(FULL_ANNOTATION_HEADER)
            short_header_str = ";".join(short_header)
            outfile_debug.write(f"{full_annotation_header_str}\n")
            outfile.write(f"{short_header_str}\n")
            for result_line in results_lists:
                # write everything to the debug files
                res_str = ";".join(map(str, result_line))
                ftid = result_line[header_dict["FTID"]]
                fusion_type = result_line[header_dict["type"]]
                ctx_id = result_line[header_dict["context_sequence_id"]]
                ctx_bp = result_line[header_dict["context_sequence_bp"]]
                ctx_wt_bp = result_line[header_dict["context_sequence_wt1_bp"]]
                ctx_wt2_bp = result_line[header_dict['context_sequence_wt2_bp']]
                outfile_debug.write(f"{res_str}\n")
                SeqIO.write(
                    SeqRecord(
                        result_line[header_dict["fusion_transcript"]],
                        id=f"{ftid}_ft_{bp_in_fusion_nt}",
                        name="",
                        description="",
                    ),
                    tfasta,
                    "fasta",
                )
                SeqIO.write(
                    SeqRecord(
                        result_line[header_dict["fusion_peptide"]],
                        id=f"{ftid}_{bp_in_fusion_aa}_{fusion_type}",
                        name="",
                        description="",
                    ),
                    pfasta,
                    "fasta",
                )

                # but filter all annotation biases (doupled ftid_plus)
                # and fusion with transcripts which have a TSL of 4 or 5
                # (keeping 1-3, 6(liftover), NA)
                if (
                    result_line[header_dict["annotation_bias"]]
                    or result_line[header_dict["wt1_TSL"]]
                    in tsl_filter_level.split(",")
                    or result_line[header_dict["wt2_TSL"]]
                    in tsl_filter_level.split(",")
                ):
                    continue
                # data for csv file
                data_str = ";".join(
                            [
                                str(result_line[header_dict[col_name]])
                                for col_name in short_header
                            ]
                )
                outfile.write(f"{data_str}\n")
                # context data for star genome generation
                count_context_seqs += 3
                summed_context_seqs_len += sum(
                    map(
                        len,
                        [
                            result_line[header_dict["context_sequence"]],
                            result_line[header_dict["context_sequence_wt1"]],
                            result_line[header_dict["context_sequence_wt2"]],
                        ],
                    )
                )
                # sequences for requantification
                SeqIO.write(
                    SeqRecord(
                        result_line[header_dict["context_sequence"]],
                        id=f"{ftid}_{ctx_id}_{ctx_bp}_ft",
                        name="",
                        description="",
                    ),
                    cfasta,
                    "fasta",
                )
                SeqIO.write(
                    SeqRecord(
                        result_line[header_dict["context_sequence_wt1"]],
                        id=f"{ftid}_{ctx_id}_{ctx_wt_bp}_wt1",
                        name="",
                        description="",
                    ),
                    cfasta,
                    "fasta",
                )
                SeqIO.write(
                    SeqRecord(
                        result_line[header_dict["context_sequence_wt2"]],
                        id=f"{ftid}_{ctx_id}_{ctx_wt2_bp}_wt2",
                        name="",
                        description="",
                    ),
                    cfasta,
                    "fasta",
                )
        # write fasta.info file for star genome generate
        with open(f"{context_seqs_file}.fasta.info", "w", encoding="utf8") as ifasta:
            ifasta.write(f"{count_context_seqs}\n")
            ifasta.write(f"{summed_context_seqs_len}\n")


def main():
    """main stub"""
    parser = ArgumentParser(description="Annotates detected fusion events")
    parser.add_argument(
        "--detected_fusions",
        dest="detected_fusions",
        help="Detected_Fusions.csv file from fusiontoolparser",
        required=True,
    )
    parser.add_argument(
        "--annotation_db",
        dest="annotation_db",
        help="gff3 derived annotation database from gffutils",
        required=True,
    )
    parser.add_argument(
        "--tsl_info",
        dest="tsl_info",
        help="File containing transcript support level (TSL). \
            Should have been generated by gtf2tsl.py",
        required=True,
    )
    parser.add_argument(
        "--genome_fasta",
        dest="genome_fasta",
        help="Genome fasta file for sequence grepping",
        required=True,
    )
    parser.add_argument(
        "--cis_near_dist",
        dest="cis_near_dist",
        type=int,
        help="Genomic distance after which a cis_near becomes a cis_far fusion",
        default=400,
    )
    parser.add_argument(
        "--context_seq_length",
        dest="context_seq_len",
        type=int,
        help="The maximal length of sequence context (exon) at either site of the bp",
        default=1000000,
    )
    parser.add_argument(
        "--tsl_filter_level",
        dest="tsl_filter_level",
        help="The TSL level to exclude from the final tables. \
            Should be comma separated w/o spaces. Possible values include 1-5 and NA",
        required=True,
    )
    parser.add_argument(
        "--out_csv",
        dest="out_csv",
        help="Path to output CSV with annotated fusions",
        required=True,
    )
    args = parser.parse_args()

    fusannot = FusionAnnotation(
        args.annotation_db, args.detected_fusions, args.tsl_info
    )
    fusannot.run(
        args.out_csv,
        int(args.cis_near_dist),
        args.genome_fasta,
        int(args.context_seq_len),
        args.tsl_filter_level,
    )

if __name__ == "__main__":
    main()
