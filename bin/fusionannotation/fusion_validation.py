"""
This module provides methods to determine fusion type.
"""

def define_type(
        cis_near_distance, bp1_chr, bp1_pos, bp1_strand, bp2_chr, bp2_pos, bp2_strand
    ):
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


def get_fusion_feature_list(
        bp1_pos,
        bp2_pos,
        bp1_strand,
        bp2_strand,
        bp1_feature_pos_list,
        bp2_feature_pos_list,
    ):
    """
    Based on the breakpoints, the strand and the complete feature positions 
    of the involved genes, return only those feature positions which will remain in the fusion.
    """
    bp1_feature_fus_list = []
    bp2_feature_fus_list = []

    # get fusion partner 1 cds
    if not len(bp1_feature_pos_list) == 0:
        if (
            bp1_strand == "+"
        ):  # on the "+" strand, we need everything LEFT of the bp for fusion gene partner 1
            bp1_feature_fus_list = [
                feature for feature in bp1_feature_pos_list if bp1_pos >= feature[0]
            ]
            if not len(bp1_feature_fus_list) == 0:
                # delete the last cds (i.e. last before the bp) if it starts with the breakpoint
                if bp1_feature_fus_list[-1][0] == bp1_pos:
                    bp1_feature_fus_list = bp1_feature_fus_list[:-1]
                else:  # or shorten the cds to the bp if it is within it
                    bp1_feature_fus_list[-1] = (
                        bp1_feature_fus_list[-1][0],
                        bp1_pos,
                    )
        else:  # on the "-" strand, we need everything RIGHT of the bp for fusion gene partner 1
            bp1_feature_fus_list = [
                feature for feature in bp1_feature_pos_list if bp1_pos <= feature[1]
            ]
            if not len(bp1_feature_fus_list) == 0:
                # delete first cds (i.e. first after the bp) if it starts with the breakpoint
                if bp1_feature_fus_list[0][1] == bp1_pos:
                    bp1_feature_fus_list = bp1_feature_fus_list[1:]
                else:
                    bp1_feature_fus_list[0] = (bp1_pos, bp1_feature_fus_list[0][1])

    # get fusion partner 2 cds
    if not len(bp2_feature_pos_list) == 0:
        if (
            bp2_strand == "+"
        ):  # on the "+" strand, we need everything RIGHT of the bp for fusion gene partner 2
            bp2_feature_fus_list = [
                feature for feature in bp2_feature_pos_list if bp2_pos <= feature[1]
            ]
            if not len(bp2_feature_fus_list) == 0:
                # delete first cds (i.e. first after the bp) if it starts with the breakpoint
                if bp2_feature_fus_list[0][1] == bp2_pos:
                    bp2_feature_fus_list = bp2_feature_fus_list[1:]
                else:
                    bp2_feature_fus_list[0] = (bp2_pos, bp2_feature_fus_list[0][1])
        else:  # on the "-" strand, we need everything LEFT of the bp for fusion gene partner 2
            bp2_feature_fus_list = [
                feature for feature in bp2_feature_pos_list if bp2_pos >= feature[0]
            ]
            if not len(bp2_feature_fus_list) == 0:
                # delete the last cds (i.e. last before the bp) if it starts with the breakpoint
                if bp2_feature_fus_list[-1][0] == bp2_pos:
                    bp2_feature_fus_list = bp2_feature_fus_list[:-1]
                else:  # or shorten the cds to the bp if it is within it
                    bp2_feature_fus_list[-1] = (
                        bp2_feature_fus_list[-1][0],
                        bp2_pos,
                    )
    return bp1_feature_fus_list, bp2_feature_fus_list
