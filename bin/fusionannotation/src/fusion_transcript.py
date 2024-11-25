"""
Module for FusionTranscript class.
"""

# pylint: disable=E0401
from .breakpoint import Breakpoint
from .transcript import Transcript


EXON_BOUNDARY_LEFT = "left_boundary"
EXON_BOUNDARY_RIGHT = "right_boundary"
EXON_BOUNDARY_WITHIN = "within"
EXON_BOUNDARY_OUTSIDE = "outside"
EXON_BOUNDARY_NA = "NA"
FUSION_TYPE_TRANS = "trans"


def get_involved_features(features: list, bp: Breakpoint, reverse: bool = False) -> list:
    """
    Based on the breakpoints, the strand and the complete feature positions 
    of the involved genes, return only those feature positions which will remain in the fusion.
    """
    bp_feature_fus_list = []
    # Does the strand relate to the exon strand or the breakpoint strand?
    # get fusion partner cds
    if not len(features) == 0:
        if (
            bp.strand == "+" and not reverse or bp.strand == "-" and reverse
        ):  # on the "+" strand, we need everything LEFT of the bp for fusion gene partner 1
            bp_feature_fus_list = [
                feature for feature in features if bp.pos >= feature.start
            ]
            if not len(bp_feature_fus_list) == 0:
                # delete the last cds (i.e. last before the bp) if it starts with the breakpoint
                if bp_feature_fus_list[-1].start == bp.pos:
                    bp_feature_fus_list = bp_feature_fus_list[:-1]
                else:  # or shorten the cds to the bp if it is within it
                    bp_feature_fus_list[-1].pos = bp.pos
        else:  # on the "-" strand, we need everything RIGHT of the bp for fusion gene partner 1
            bp_feature_fus_list = [
                feature for feature in features if bp.pos <= feature.stop
            ]
            if not len(bp_feature_fus_list) == 0:
                # delete first cds (i.e. first after the bp) if it starts with the breakpoint
                if bp_feature_fus_list[0].stop == bp.pos:
                    bp_feature_fus_list = bp_feature_fus_list[1:]
                else:
                    bp_feature_fus_list[0].pos = bp.pos

    return bp_feature_fus_list


class FusionTranscript:
    """
    Class to store fusion transcript information
    """
    def __init__(
        self,
        transcript_1: Transcript,
        transcript_2: Transcript,
        bp1: Breakpoint,
        bp2: Breakpoint,
        cis_near_distance: int = 1000000
    ):
        self.transcript_1 = transcript_1
        self.transcript_2 = transcript_2
        self.bp1 = bp1
        self.bp2 = bp2
        self.frame = self.get_fusion_frame()
        # Why is the reverse flag set to a fixed value? --> Should this represent the orientation of the fusion? (LR, LL, RL, RR)
        self.exons_transcript_1 = get_involved_features(self.transcript_1.exons, self.bp1, reverse=False)
        self.exons_transcript_2 = get_involved_features(self.transcript_2.exons, self.bp2, reverse=True)
        self.cds_transcript_1 = get_involved_features(self.transcript_1.cds, self.bp1, reverse=False)
        self.cds_transcript_2 = get_involved_features(self.transcript_2.cds, self.bp2, reverse=True)
        self.fusion_type = self.determine_fusion_type(cis_near_distance)


    def __repr__(self):
        return f"FusionTranscript(id={repr(self.get_ftid())}, " \
               f"bpid={repr(self.get_bpid())})"


    def get_bpid(self):
        """Return the breakpoint ID"""
        return f"{self.bp1}_{self.bp2}"


    def get_ftid(self):
        """Return the fusion transcript ID"""
        return "_".join(
            [
                self.transcript_1.gene_name,
                str(self.bp1),
                self.transcript_1.transcript_id.strip("transcript:"),
                self.transcript_2.gene_name,
                str(self.bp2),
                self.transcript_2.transcript_id.strip("transcript:"),
            ]
        )


    def get_fusion_gene_name(self):
        """Return the fusion gene name"""
        return "_".join(
            [
                self.transcript_1.gene_name,
                self.transcript_2.gene_name,
            ]
        )

    # Is this used after all suspect transcripts have been annotated?
    def set_flags(self, suspect_transcripts: dict):
        """Set flags for suspect transcripts"""
        self.transcript_1.flags = suspect_transcripts[self.transcript_1.transcript_id]
        self.transcript_2.flags = suspect_transcripts[self.transcript_2.transcript_id]


    def annotate_same_chrom(self, cis_near_distance: int) -> str:
        """Annotates the fusion type for same chromosome and same strand fusions.

        Args:
            cis_near_distance (int): Distance to classify as cis_near fusions

        Returns:
            str: Type for same chromosome and same strand fusions
        """
        distance = 0
        if self.bp1.strand == "+":
            distance = self.bp2.pos - self.bp1.pos
        else:
            distance = self.bp1.pos - self.bp2.pos
        if distance < 0:
            return "cis_trans"
        if 0 < distance <= cis_near_distance:
            return "cis_near"
        return "cis_far"


    def determine_fusion_type(self, cis_near_distance: int) -> str:
        """Determine the fusion type based on the location and orientation of the fusion partners"""
        # Fusion type:
        # -	"trans" (bp of fusion partners on different chromosomes, but with same strand)
        # -	"trans_inv" (bp of fusion partners on different chromosomes and with different strand)
        # -	"cis_inv" (bp of fusion partners on same chromosome, but with different strand)
        # -	"cis_trans" (bp of fusion partners on same chromosome and with same strand,
        #       but 2nd fusion partner is before (strandwise) 1st fusion partner)
        # -	"cis_far" (bp of fusion partners on same chromosome, same strand,
        #       1st before 2nd, but more than 1Mb apart from each other (genomic distance))
        # -	"cis_near" (bp of fusion partners on same chromosome, same strand,
        #       1st before 2nd and closer than 1Mb apart from each other (genomic distance))
        same_chrom = self.bp1.chrom == self.bp2.chrom
        same_strand = self.bp1.strand == self.bp2.strand

        fusion_types = {
            (True, True): self.annotate_same_chrom(cis_near_distance),
            (True, False): "cis_inv",
            (False, True): "trans",
            (False, False): "trans_inv"
        }

        return fusion_types.get((same_chrom, same_strand), "no_match")


    def get_combined_boundary(self) -> str:
        """Check whether boundaries of both fusion partners are in the correct orientation"""

        is_on_boundary = any((EXON_BOUNDARY_LEFT, EXON_BOUNDARY_RIGHT))
        stranded = any(("+", "-"))

        boundary_dict = {
            ("+", "+", EXON_BOUNDARY_RIGHT, EXON_BOUNDARY_LEFT): "both",
            ("+", "-", EXON_BOUNDARY_RIGHT, EXON_BOUNDARY_RIGHT): "both",
            ("-", "+", EXON_BOUNDARY_LEFT, EXON_BOUNDARY_LEFT): "both",
            ("-", "-", EXON_BOUNDARY_LEFT, EXON_BOUNDARY_RIGHT): "both",
            (stranded, stranded, is_on_boundary, not is_on_boundary): "5prime",
            (stranded, stranded, not is_on_boundary, is_on_boundary): "3prime"
        }
        return boundary_dict.get(
            (self.bp1.strand, self.bp2.strand, self.bp1.exon_boundary, self.bp2.exon_boundary),
            "no_match"
        )


    def get_fusion_frame(self) -> str:
        """Combine frame info from bp1 and bp2"""

        if self.transcript_1.frame_at_bp == -1:
            return "no_frame"
        if self.transcript_2.frame_at_bp == -1:
            return "neo_frame"
        if self.transcript_1.frame_at_bp == self.transcript_2.frame_at_bp:
            return "in_frame"
        return "out_frame"


    def get_exon_nr(self) -> int:
        """Get the number of exons in the fusion transcript"""
        if self.frame == "in_frame":
            return (len(self.cds_transcript_1) +
                    len(self.cds_transcript_2))
        return (len(self.cds_transcript_1) +
                len(self.exons_transcript_2))


    def has_overlapping_transcripts(self) -> bool:
        """Check if exons of wt1 and wt2 overlap.
        Since both exon lists are sorted, we can simply compare the first 
        and last exon of each list.

        Args:
            wt1_exons (list): _description_
            wt2_exons (list): _description_
            fusion_type (str): _description_

        Returns:
            bool: True if there is an overlap, False otherwise
        """
        # Skip check for translocations
        if self.fusion_type.startswith(FUSION_TYPE_TRANS):
            return False
        # If there are no exons in either of the transcripts, there is no overlap
        if len(self.transcript_1.exons) == 0 or len(self.transcript_2.exons) == 0:
            return False
        # end of wt1 locus must be < start of wt2 locus or start of wt1 locus > end of wt2 locus
        if (
            self.transcript_1.exons[-1].stop < self.transcript_2.exons[0].start
            or self.transcript_1.exons[0].start > self.transcript_2.exons[-1].stop
        ):
            return False
        return True
