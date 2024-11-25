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
@author: TRON (PASO)
@version: 241118
"""

from argparse import ArgumentParser
import logging

# pylint: disable=E0401
from fusionannotation.src.breakpoint import Breakpoint
from fusionannotation.src.io_methods import load_detected_fusions
from fusionannotation.src.io_methods import load_tsl_data
from fusionannotation.src.io_methods import load_genomic_data
from fusionannotation.src.feature_validation import get_exon_cds_overlap
from fusionannotation.src.fusion_transcript import FusionTranscript
from fusionannotation.src.gff_db_controller import GffDbController
from fusionannotation.src.output_handler import OutputHandler
from fusionannotation.src.result_handler import ResultHandler
from fusionannotation.src.transcript import Transcript


FORMAT = '%(asctime)s - %(message)s'
logging.basicConfig(format=FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)


class FusionAnnotator:
    """Annotation of predicted fusion genes soley based on the breakpoint information"""

    def __init__(self, db_file, bp_file, tsl_file):
        """Parameter initialization"""
        logger.info("Loading fusion breakpoints into dict")
        self.bp_dict = load_detected_fusions(bp_file)
        logger.info("Loading GFF data into database object")
        self.db = GffDbController(db_file)
        logger.info("Loading TSL data into dict")
        self.tsl_dict = load_tsl_data(tsl_file)
        self.suspect_transcripts = {}


    def annotate_bp(self, bp: Breakpoint, exon_id: str) -> Transcript:
        """
        Annotate breakpoint of a pair
        """
        # get id and additional infos from exon parents transcript and gene
        transcript = self.db.get_parent_transcript(exon_id)
        # create a set to collect "suspicious findings" for transcripts
        if transcript.transcript_id not in self.suspect_transcripts:
            self.suspect_transcripts[transcript.transcript_id] = set()
        # get the complete set of exons, corresponding cds
        # (cds at the same index are NA if there is none in
        # the exon) and start frames of the cds
        exons = self.db.get_features_from_transcript(transcript.transcript_id, "exon")
        cds = self.db.get_features_from_transcript(transcript.transcript_id, "CDS")
        (cds_pos_list, trans_flags) = get_exon_cds_overlap(
            exons, cds
        )
        self.suspect_transcripts[transcript.transcript_id].update(trans_flags)
        transcript.set_exons(exons)
        # get the frame at the start of the first cds and at the breakpoint
        frame_at_start, frame_at_bp = bp.get_frame(
            cds_pos_list
        )
        transcript.set_frame_at_start(frame_at_start)
        transcript.set_frame_at_bp(frame_at_bp)
        if frame_at_start > 0:
            # flag the transcript if its frame at the start of the first cds is > 0
            self.suspect_transcripts[transcript.transcript_id].add("cds start > 0")
        # remove "NA's" from the cds list
        cds_pos_list = [x for x in cds_pos_list if x]
        transcript.set_cds(cds_pos_list)
        if len(cds_pos_list) == 1:
            # flag the transcript if it contains a single cds
            self.suspect_transcripts[transcript.transcript_id].add("Only 1 cds in bp1")
        # add pre-loaded tsl info
        # The "correct" tsl for a missing trans_id would be "NA",
        # I'm using "6" in order to allow differentiation from
        # existing, "normal" tsl's which range from 1-5 and are NA
        # for transcripts of unknown support level
        tsl = self.tsl_dict.get(transcript.transcript_id, "6")
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

        if not bpid == "12:11869969:+_15:87940753:-":
            return []
        # get features overlapping with the bp from the db
        (bp1_exons, bp1_cds) = self.db.get_exons_cds_overlapping_position(bp1)
        (bp2_exons, bp2_cds) = self.db.get_exons_cds_overlapping_position(bp2)
        # consider all possible combination between transcript1 to transcript2 fusions
        for bp1_efeature, bp1_cfeature in zip(bp1_exons, bp1_cds):
            bp1.exon_boundary = bp1.get_boundary(bp1_efeature)
            bp1.cds_boundary = bp1.get_boundary(bp1_cfeature)
            wt1 = self.annotate_bp(bp1, bp1_efeature.id)
            # do the same stuff for transcript 2 and (if applicable)
            # combine the information from both transcripts
            for bp2_efeature, bp2_cfeature in zip(bp2_exons, bp2_cds):
                bp2.exon_boundary = bp2.get_boundary(bp2_efeature)
                bp2.cds_boundary = bp2.get_boundary(bp2_cfeature)
                wt2 = self.annotate_bp(bp2, bp2_efeature.id)
                fusion_transcript = FusionTranscript(wt1, wt2, bp1, bp2, cis_near_distance)
                # check whether loci of wt1/wt2 overlap
                # Why flag only wt1? If wt1 and wt2 overlap, then wt2 also overlaps wt1
                # It makes more sense to flag fusion transcript instead of the wildtype transcript,
                # since the wildtype transcript can be used in multiple fusion transcripts
                if fusion_transcript.has_overlapping_transcripts():
                    self.suspect_transcripts[wt1.transcript_id].add(
                        "wt1/wt2 exon overlap"
                    )
                fusion_transcript.set_flags(self.suspect_transcripts)
                result_list.append(fusion_transcript)
        
        # An improvement may be to add flags after the annotation process to catch all suspect_transcripts.
        # for fusion_transcript in result_list:
        #     fusion_transcript.set_flags(self.suspect_transcripts)

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

        logger.info("Running annotation pipeline")

        fusion_transcripts = []

        for bpid in sorted(self.bp_dict):
            fusion_transcripts_bpid = self.annotate_bpid(bpid, cis_near_distance)
            fusion_transcripts.extend(fusion_transcripts_bpid)

        logger.info("Loading genomic data from %s", genome_fasta)
        cds_seqs = load_genomic_data(genome_fasta, fusion_transcripts)

        result_handler = ResultHandler(cds_seqs, context_seq_len)
        result_list = result_handler.generate_result_list(
            fusion_transcripts
        )
        output_handler = OutputHandler(result_list, context_seqs_file, tsl_filter_level)
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

    fusannot = FusionAnnotator(
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
