"""
Module for FusionTranscript class.
"""
import copy
from typing import Union, List
# pylint: disable=E0401
from .breakpoint import Breakpoint
from .transcript import Transcript
from .exon import Exon
from .cds import CDS


EXON_BOUNDARY_LEFT = "left_boundary"
EXON_BOUNDARY_RIGHT = "right_boundary"
EXON_BOUNDARY_WITHIN = "within"
EXON_BOUNDARY_OUTSIDE = "outside"
EXON_BOUNDARY_NA = "NA"
FUSION_TYPE_TRANS = "trans"


def get_involved_features(
    features: List[Union[CDS, Exon]], bp: Breakpoint, reverse: bool = False
) -> List[Union[CDS, Exon]]:
    """
    Return truncated feature list (prefix or suffix) relative to breakpoint.
    Does not mutate the original feature objects.
    reverse:
      False -> keep 5' (prefix) part for '+' strand (left), 3' part for '-' strand
      True  -> inverted (used for partner 2 in current logic)
    """
    if not features:
        return []
    feats = copy.deepcopy(features)

    keep = []
    if (bp.strand == "+" and not reverse) or (bp.strand == "-" and reverse):
        # keep LEFT of breakpoint
        keep = [f for f in feats if bp.pos >= f.start]
        if keep:
            if keep[-1].start == bp.pos:
                keep = keep[:-1]
            elif keep[-1].start < bp.pos <= keep[-1].stop:
                keep[-1].stop = bp.pos
    else:
        # keep RIGHT of breakpoint
        keep = [f for f in feats if bp.pos <= f.stop]
        if keep:
            if keep[0].stop == bp.pos:
                keep = keep[1:]
            elif keep[0].start <= bp.pos < keep[0].stop:
                keep[0].start = bp.pos
    return keep


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
        self.transcript_1 = copy.deepcopy(transcript_1)
        self.transcript_2 = copy.deepcopy(transcript_2)
        self.bp1 = copy.deepcopy(bp1)
        self.bp2 = copy.deepcopy(bp2)

        # Truncated feature sets (do not overwrite originals)
        self.exons_transcript_1 = get_involved_features(self.transcript_1.exons, self.bp1, reverse=False)
        self.exons_transcript_2 = get_involved_features(self.transcript_2.exons, self.bp2, reverse=True)
        self.cds_transcript_1 = get_involved_features(self.transcript_1.cds, self.bp1, reverse=False)
        self.cds_transcript_2 = get_involved_features(self.transcript_2.cds, self.bp2, reverse=True)

        self.frame = self.get_fusion_frame()
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
        """Determine the combined boundary type of both fusion partners."""

        # Check if either breakpoint is on a valid exon boundary
        is_on_boundary_1 = self.bp1.exon_boundary in [EXON_BOUNDARY_LEFT, EXON_BOUNDARY_RIGHT]
        is_on_boundary_2 = self.bp2.exon_boundary in [EXON_BOUNDARY_LEFT, EXON_BOUNDARY_RIGHT]

        # Define valid boundary combinations for "both"
        valid_combinations = {
            ("+", "+"): (EXON_BOUNDARY_RIGHT, EXON_BOUNDARY_LEFT),
            ("+", "-"): (EXON_BOUNDARY_RIGHT, EXON_BOUNDARY_RIGHT),
            ("-", "+"): (EXON_BOUNDARY_LEFT, EXON_BOUNDARY_LEFT),
            ("-", "-"): (EXON_BOUNDARY_LEFT, EXON_BOUNDARY_RIGHT),
        }

        # Check if the current combination matches any valid "both" boundary
        if (self.bp1.strand, self.bp2.strand) in valid_combinations:
            expected_boundaries = valid_combinations[(self.bp1.strand, self.bp2.strand)]
            if (self.bp1.exon_boundary, self.bp2.exon_boundary) == expected_boundaries:
                return "both"

        # Determine if the fusion is 5' or 3' based on boundary presence
        if is_on_boundary_1 and not is_on_boundary_2:
            return "5prime"
        if not is_on_boundary_1 and is_on_boundary_2:
            return "3prime"

        # Default case if no match is found
        return "no_match"


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
        """Number of contributing (truncated) exons/CDS."""
        if self.frame == "in_frame":
            return len(self.cds_transcript_1) + len(self.cds_transcript_2)
        return len(self.cds_transcript_1) + len(self.exons_transcript_2)


    def has_overlapping_transcripts(self) -> bool:
        """Overlap check on original (untruncated) loci."""
        if self.fusion_type.startswith(FUSION_TYPE_TRANS):
            return False
        if not self.transcript_1.exons or not self.transcript_2.exons:
            return False
        if (self.transcript_1.exons[-1].stop < self.transcript_2.exons[0].start or
            self.transcript_1.exons[0].start > self.transcript_2.exons[-1].stop):
            return False
        return True
