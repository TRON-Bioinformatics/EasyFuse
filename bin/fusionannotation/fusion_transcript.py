"""
Module for FusionTranscript class.
"""

# pylint: disable=E0401
from breakpoint import Breakpoint
from transcript import Transcript


def get_involved_exons(features: list, bp: Breakpoint, reverse: bool = False) -> list:
    """
    Based on the breakpoints, the strand and the complete feature positions 
    of the involved genes, return only those feature positions which will remain in the fusion.
    """
    bp_feature_fus_list = []

    # get fusion partner 1 cds
    if not len(features) == 0:
        if (
            bp.strand == "+" and not reverse
        ):  # on the "+" strand, we need everything LEFT of the bp for fusion gene partner 1
            bp_feature_fus_list = [
                feature for feature in features if bp.pos >= feature[0]
            ]
            if not len(bp_feature_fus_list) == 0:
                # delete the last cds (i.e. last before the bp) if it starts with the breakpoint
                if bp_feature_fus_list[-1][0] == bp.pos:
                    bp_feature_fus_list = bp_feature_fus_list[:-1]
                else:  # or shorten the cds to the bp if it is within it
                    bp_feature_fus_list[-1] = (
                        bp_feature_fus_list[-1][0],
                        bp.pos,
                    )
        else:  # on the "-" strand, we need everything RIGHT of the bp for fusion gene partner 1
            bp_feature_fus_list = [
                feature for feature in features if bp.pos <= feature[1]
            ]
            if not len(bp_feature_fus_list) == 0:
                # delete first cds (i.e. first after the bp) if it starts with the breakpoint
                if bp_feature_fus_list[0][1] == bp.pos:
                    bp_feature_fus_list = bp_feature_fus_list[1:]
                else:
                    bp_feature_fus_list[0] = (bp.pos, bp_feature_fus_list[0][1])

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
        bp2: Breakpoint
    ):
        self.transcript_1 = transcript_1
        self.transcript_2 = transcript_2
        self.bp1 = bp1
        self.bp2 = bp2
        self.frame = self.get_combined_frame()
        self.exons_transcript_1 = self.get_fusion_exons_transcript_1()
        self.exons_transcript_2 = self.get_fusion_exons_transcript_2()
        self.cds_transcript_1 = self.get_fusion_cds_transcript_1()
        self.cds_transcript_2 = self.get_fusion_cds_transcript_2()


    def get_ftid(self):
        """Return the fusion transcript ID"""
        return "_".join(
            [
                self.transcript_1.gene_name,
                self.bp1,
                self.transcript_1.transcript_id.split(":")[1],
                self.transcript_2.gene_name,
                self.bp2,
                self.transcript_2.transcript_id.split(":")[1],
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


    def get_fusion_type(self, cis_near_distance: int) -> str:
        """Define the fusion type based on the location and orientation of the fusion partners"""
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

        if self.bp1.chrom == self.bp2.chrom:
            if self.bp1.strand == self.bp2.strand:
                if self.bp1.strand == "+":
                    if self.bp2.pos - self.bp1.pos < 0:
                        return "cis_trans"
                    if self.bp2.pos - self.bp1.pos > cis_near_distance:
                        return "cis_far"
                    return "cis_near"
                else:
                    if self.bp1.pos - self.bp2.pos < 0:
                        return "cis_trans"
                    if self.bp1.pos - self.bp2.pos > cis_near_distance:
                        return "cis_far"
                    return "cis_near"
            return "cis_inv"
        else:
            if self.bp1.strand == self.bp2.strand:
                return "trans"
            return "trans_inv"


    def get_combined_frame(self) -> str:
        """Combine frame info from bp1 and 2"""

        if self.transcript_1.frame == -1:
            return "no_frame"
        if self.transcript_2.frame == -1:
            return "neo_frame"
        if self.transcript_1.frame == self.transcript_2.frame:
            return "in_frame"
        return "out_frame"


    def get_fusion_cds_transcript_1(self):
        """Return the fusion CDS for transcript 1"""
        return get_involved_exons(self.transcript_1.cds_pos_list, self.bp1)


    def get_fusion_cds_transcript_2(self):
        """Return the fusion CDS for transcript 2"""
        return get_involved_exons(self.transcript_2.cds_pos_list, self.bp2, reverse=True)


    def get_fusion_exons_transcript_1(self):
        """Return the fusion exons for transcript 1"""
        return get_involved_exons(self.transcript_1.exon_pos_list, self.bp1)


    def get_fusion_exons_transcript_2(self):
        """Return the fusion exons for transcript 2"""
        return get_involved_exons(self.transcript_2.exon_pos_list, self.bp2, reverse=True)


    def get_exon_nr(self) -> int:
        """Get the number of exons in the fusion transcript"""
        if self.frame == "in_frame":
            return (len(self.cds_transcript_1) +
                    len(self.cds_transcript_2))
        return (len(self.cds_transcript_1) +
                len(self.exons_transcript_2))
