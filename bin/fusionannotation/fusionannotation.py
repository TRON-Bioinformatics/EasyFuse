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
import gffutils # type: ignore
from logzero import logger # type: ignore


from breakpoint import Breakpoint

from io_methods import load_detected_fusions
from io_methods import load_tsl_data
from boundary_validation import get_combined_boundary
from exon_validation import check_exon_overlap
from exon_validation import get_wt_codings
from exon_validation import get_parents
from frame_validation import get_frame
from fusion_transcript import FusionTranscript
from output_handler import OutputHandler
from result_handler import update_results




class FusionAnnotation:
    """Annotation of predicted fusion genes soley based on the breakpoint information"""

    def __init__(self, db_file, bp_file, tsl_file):
        """Parameter initialization"""
        self.bp_dict = load_detected_fusions(bp_file)
        self.db = gffutils.FeatureDB(db_file)
        self.tsl_dict = load_tsl_data(tsl_file)
        self.cds_seq_dict = {}
        self.suspect_transcripts = {}


    def get_tsl(self, trans_id: str) -> str:
        """Get the transcript support level (tsl) from the pre-loaded dict"""
        # The "correct" tsl for a missing trans_id would be "NA",
        # I'm using "6" in order to allow differentiation from
        # existing, "normal" tsl's which range from 1-5 and are NA
        # for transcripts of unknown support
        if not trans_id in self.tsl_dict:
            return "6"
        return self.tsl_dict[trans_id]


    def fill_seq_lookup_dict(self, chrom: str, cds_list: list):
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


    def fill_seq_lookup_dict_fusion_transcript(self, ft: FusionTranscript):
        """Collect positions for which we need to grep sequences."""
        self.fill_seq_lookup_dict(ft.bp1.chrom, ft.transcript_1.exon_pos_list)
        self.fill_seq_lookup_dict(ft.bp2.chrom, ft.transcript_2.exon_pos_list)
        self.fill_seq_lookup_dict(ft.bp1.chrom, ft.exons_transcript_1)
        self.fill_seq_lookup_dict(ft.bp2.chrom, ft.exons_transcript_2)
        self.fill_seq_lookup_dict(ft.bp1.chrom, ft.transcript_1.cds_pos_list)
        self.fill_seq_lookup_dict(ft.bp2.chrom, ft.transcript_2.cds_pos_list)
        self.fill_seq_lookup_dict(ft.bp1.chrom, ft.cds_transcript_1)
        self.fill_seq_lookup_dict(ft.bp2.chrom, ft.cds_transcript_2)


    def grep_candidate_seqs(self, genome_fasta: str):
        """
        Grep chromosomes and extract cds sequences from the @cds_seq_dict
        urla note: the max(..) check in the slicing is probably superfluous
        as CDS annotations should not start at chromosome position 0
        (which would become -1 in @get_frame()). Nevertheless, negative values would
        completely mess up the slicing...
        An additinal check at the end is not required as
        slicing stops anyway at the end of sequence
        """

        for record in SeqIO.parse(genome_fasta, "fasta"):
            if record.id in self.cds_seq_dict:
                logger.info(
                    "Grepping candidate regions from chromosome %s", record.id
                )

                self.cds_seq_dict[record.id][1] = [
                    record.seq[max(0, val[0] - 1) : val[1]]
                    for val in self.cds_seq_dict[record.id][0]
                ]


    def annotate_bp(self, bp: Breakpoint, bp_efeature: object):
        """
        Annotate the first breakpoint of a pair
        """
        # get id and additional infos from exon parents transcript and gene
        transcript = get_parents(self.db, bp_efeature)
        # create a set to collect "suspicious findings" for transcripts
        if transcript.transcript_id not in self.suspect_transcripts:
            self.suspect_transcripts[transcript.transcript_id] = set()
        # get the complete set of exons, corresponding cds
        # (cds at the same index are NA if there is none in
        # the exon) and start frames of the cds
        (exon_pos_list, cds_pos_list, cds_frame_list) = get_wt_codings(
            self.db,
            self.suspect_transcripts,
            transcript.transcript_id
        )
        transcript.set_exon_pos_list(exon_pos_list)
        # get the frame at the start of the first cds and at the breakpoint
        frame_at_start, frame_at_bp = get_frame(
            bp, cds_pos_list, cds_frame_list
        )
        transcript.set_frame_at_bp(frame_at_bp)
        if frame_at_start > 0:
            # flag the transcript if its frame at the start of the first cds is > 0
            self.suspect_transcripts[transcript.transcript_id].add("cds start > 0")
        # remove "NA's" from the cds list
        cds_pos_list = [x for x in cds_pos_list if not x == "NA"]
        transcript.set_cds_pos_list(cds_pos_list)
        if len(cds_pos_list) == 1:
            # flag the transcript if it contains a single cds
            self.suspect_transcripts[transcript.transcript_id].add("Only 1 cds in bp1")
        # add pre-loaded tsl info
        tsl = self.get_tsl(transcript.transcript_id.split(":")[1])
        transcript.set_tsl(tsl)
        return transcript


    def annotate_bpid(self, bpid: str, cis_near_distance: int):
        """
        Annotate a single breakpoint pair
        """
        logger.debug("Annotating '%s'...", bpid)
        result_list = []
        # split breakpoint info
        bp1 = self.bp_dict[bpid][0]
        bp2 = self.bp_dict[bpid][1]

        # query db to get exons and matched cds (if available) overlapping
        # the bp as gffutils feature objects
        bp1_exons, bp1_cds = bp1.get_overlapping_features(self.db)
        bp2_exons, bp2_cds = bp2.get_overlapping_features(self.db)
        # consider all possible combination between transcript1 to transcript2 fusions
        for bp1_efeature, bp1_cfeature in zip(bp1_exons, bp1_cds):
            # check whether breakpoint1 is located on an exon boundary
            exon_boundary_1 = bp1.get_boundary(bp1_efeature)
            # check whether breakpoint1 is located on a CDS boundary
            cds_boundary_1 = bp1.get_boundary(bp1_cfeature)
            wt1 = self.annotate_bp(bp1, bp1_efeature)
            # do the same stuff for transcript 2 and (if applicable)
            # combine the information from both transcripts
            for bp2_efeature, bp2_cfeature in zip(bp2_exons, bp2_cds):
                # check whether breakpoint1 is located on an exon boundary
                exon_boundary_2 = bp2.get_boundary(bp2_efeature)
                # check whether breakpoint1 is located on a CDS boundary
                cds_boundary_2 = bp2.get_boundary(bp2_cfeature)
                wt2 = self.annotate_bp(bp2, bp2_efeature)
                # combined exon boundary information
                exon_boundary = get_combined_boundary(
                    exon_boundary_1, exon_boundary_2, bp1.strand, bp2.strand
                )
                fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2)
                fusion_type = fusion_transcript.get_fusion_type(cis_near_distance)

                # check whether loci of wt1/wt2 overlap
                if check_exon_overlap(
                    wt1.exon_pos_list, wt2.exon_pos_list, fusion_type
                ):
                    self.suspect_transcripts[wt1.transcript_id].add("wt1/wt2 exon overlap")

                self.fill_seq_lookup_dict_fusion_transcript(fusion_transcript)

                wt1_is_good_transcript = self.suspect_transcripts[wt1.transcript_id]
                wt2_is_good_transcript = self.suspect_transcripts[wt2.transcript_id]
                result_list.append([
                    bpid,
                    fusion_transcript.get_fusion_gene_name(),
                    bp1,
                    bp2,
                    fusion_transcript.get_ftid(),
                    "context_sequence_id",
                    "context_sequence_100_id",
                    fusion_type,
                    str(fusion_transcript.get_exon_nr()),
                    str(len(fusion_transcript.exons_transcript_1)),
                    str(len(wt1.exon_pos_list) - len(fusion_transcript.exons_transcript_2)),
                    "exon_starts",
                    "exon_ends",
                    exon_boundary_1,
                    exon_boundary_2,
                    exon_boundary,
                    str(wt1.frame),
                    str(wt2.frame),
                    fusion_transcript.frame,
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
                    bp1.chrom,
                    bp1.pos,
                    bp1.strand,
                    bp2.chrom,
                    bp2.pos,
                    bp2.strand,
                    cds_boundary_1,
                    cds_boundary_2,
                    wt1.exon_pos_list,
                    wt2.exon_pos_list,
                    fusion_transcript.exons_transcript_1,
                    fusion_transcript.exons_transcript_2,
                    wt1.cds_pos_list,
                    wt2.cds_pos_list,
                    fusion_transcript.cds_transcript_1,
                    fusion_transcript.cds_transcript_2,
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
                    wt1.transcript_biotype,
                    wt2.transcript_biotype,
                    wt1.gene_biotype,
                    wt2.gene_biotype,
                    wt1.description,
                    wt2.description,
                    wt1.frame_at_start,
                    wt2.frame_at_start,
                    wt1.tsl,
                    wt2.tsl,
                    len(wt1.exon_pos_list),
                    len(wt2.exon_pos_list),
                    len(fusion_transcript.exons_transcript_1),
                    len(fusion_transcript.exons_transcript_2),
                    len(wt1.cds_pos_list),
                    len(wt2.cds_pos_list),
                    len(fusion_transcript.cds_transcript_1),
                    len(fusion_transcript.cds_transcript_2),
                    f"{bp1.chrom}:{wt1.exon_pos_list[0][0]}:{wt1.exon_pos_list[-1][1]}",
                    f"{bp2.chrom}:{wt2.exon_pos_list[0][0]}:{wt2.exon_pos_list[-1][1]}",
                    "annotation_bias",
                ])
        return result_list


    def run(
        self,
        context_seqs_file: str,
        cis_near_distance: int,
        genome_fasta: str,
        context_seq_len: int,
        tsl_filter_level: str,
    ):
        """Run the annotation pipeline"""
        # for performance reasons (mainly seq grepping), results cannot be written
        # on the fly and everything will therefore be stored in a list of lists
        logger.info("Running annotation pipeline")


        # information in the header is complete and useful for debugging,
        # information in the short_header is identical
        # to the format of the previous
        # version of the fusion annotation and therefore used for downstream processing

        results_lists = []

        # get features overlapping with the bp from the db
        for bpid in sorted(self.bp_dict):
            result_list = self.annotate_bpid(bpid, cis_near_distance)
            results_lists.extend(result_list)

        self.grep_candidate_seqs(genome_fasta)
        results_lists = update_results(
            results_lists,
            self.cds_seq_dict,
            context_seq_len
        )
        output_handler = OutputHandler(results_lists, context_seqs_file, tsl_filter_level)
        output_handler.write_full_table()
        output_handler.write_filtered_table()
        output_handler.write_fasta()
        output_handler.write_fasta_info()
        output_handler.write_transcript_fasta()
        output_handler.write_peptide_fasta()
        logger.info("Finished annotation pipeline")


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
            Should be comma separated w/o spaces. Possible values include 1-5 and NA, \
            e.g. 1,2,NA",
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
        args.cis_near_dist,
        args.genome_fasta,
        args.context_seq_len,
        args.tsl_filter_level,
    )

if __name__ == "__main__":
    main()
