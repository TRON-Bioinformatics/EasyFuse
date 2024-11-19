"""
This module performs validation of exons.
"""

FUSION_TYPE_TRANS = "trans"

def check_exon_overlap(wt1_exon_pos_list: list, wt2_exon_pos_list: list, fusion_type: str) -> bool:
    """Returns true if exons overlap and false otherwise"""
    # exon position lists are always sorted
    # and we therefore don't need any information on the strand
    # gene on different chromosomes
    if fusion_type.startswith(FUSION_TYPE_TRANS):
        return False
    # end of wt1 locus must be < start of wt2 locus or start of wt1 locus > end of wt2 locus
    if (
        wt1_exon_pos_list[-1].stop < wt2_exon_pos_list[0].start
        or wt1_exon_pos_list[0].start > wt2_exon_pos_list[-1].stop
    ):
        return False
    return True


def get_exon_cds_overlap(exons: list, cds_list: list) -> tuple:
    """Get matched exons and cds region per transcript id"""

    cds_filt = []
    trans_flags_dict = set()
    for exon in exons:
        exon_has_cds = False
        for i, cds in enumerate(cds_list, 0):
            if cds.start >= exon.start and cds.stop <= exon.stop:
                cds_filt.append(cds)
                if cds.start != exon.start and cds.stop != exon.stop:
                    trans_flags_dict.add(f"cds != exon at no {i}")
                exon_has_cds = True
        if not exon_has_cds:
            cds_filt.append(None)

    return (cds_filt, trans_flags_dict)


def filter_cds_by_exons(exons: list, cds: list) -> list:
    """Correct positions in such a way that the same position
    in the exon/CDS list represents the same parental transcript
    exons are the scaffold as there can be exons w/o cds but not vice versa

    Args:
        exons (list): _description_
        cds (list): _description_

    Returns:
        list: _description_
    """
    cds_transcripts = [cds.transcript_id for cds in cds]
    bp_cds_new = []
    for exon in exons:
        if exon.transcript_id in cds_transcripts:
            bp_cds_new.append(cds[cds_transcripts.index(exon.transcript_id)])
        else:
            bp_cds_new.append(None)
    return bp_cds_new
