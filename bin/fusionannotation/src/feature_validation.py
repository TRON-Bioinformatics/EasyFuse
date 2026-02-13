"""
This module performs validation of exons.
"""


def get_exon_cds_overlap(exons: list, cds_list: list) -> tuple:
    """Get the CDS regions that overlap with the exons.
    Args:
        exons (list): A list of exon objects, each with 'start' and 'stop' attributes.
        cds_list (list): A list of CDS objects, each with 'start' and 'stop' attributes.

    Returns:
        tuple: A tuple containing:
            - list: A list of CDS regions that overlap with the exons.
            - set: A set of transcript flags indicating any discrepancies between the CDS and exon boundaries.
    """
    cds_filt = []
    trans_flags_dict = set()
    for exon in exons:
        exon_has_cds = False
        for i, cds in enumerate(cds_list, 0):
            if cds.start >= exon.start and cds.stop <= exon.stop:
                cds_filt.append(cds)
                if not (cds.start == exon.start and cds.stop == exon.stop):
                    trans_flags_dict.add(f"cds != exon at no {i}")
                exon_has_cds = True
            # else if cds is not in the exon --> Catch Error 
        if not exon_has_cds:
            cds_filt.append(None)

    return (cds_filt, trans_flags_dict)


def filter_cds_by_exons(exons: list, cds: list) -> list:
    """Correct positions in such a way that the same position
    in the exon/CDS list represents the same parental transcript
    exons are the scaffold as there can be exons w/o cds but not vice versa
    Make sure that CDS and exons are taken from the same transcript.

    Args:
        exons (list): A list of exon objects, each with 'start', 'stop', and 'transcript_id' attributes.
        cds (list): A list of CDS objects, each with 'start', 'stop', and 'transcript_id' attributes.

    Returns:
        list: A list of CDS objects aligned with the exons based on their transcript IDs.
    """
    cds_transcripts = [cds.transcript_id for cds in cds]
    bp_cds_new = []
    for exon in exons:
        if exon.transcript_id in cds_transcripts:
            bp_cds_new.append(cds[cds_transcripts.index(exon.transcript_id)])
        else:
            bp_cds_new.append(None)
    return bp_cds_new
