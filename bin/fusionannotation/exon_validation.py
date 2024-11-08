"""
This module performs validation of exons.
"""

def check_exon_overlap(wt1_exon_pos_list: list, wt2_exon_pos_list: list, fusion_type: str) -> bool:
    """Returns true if exons overlap and false otherwise"""
    # exon position lists are always sorted
    # and we therefore don't need any information on the strand
    # gene on different chromosomes
    if fusion_type.startswith("trans"):
        return False
    # end of wt1 locus must be < start of wt2 locus or start of wt1 locus > end of wt2 locus
    if (
        wt1_exon_pos_list[-1][1] < wt2_exon_pos_list[0][0]
        or wt1_exon_pos_list[0][0] > wt2_exon_pos_list[-1][1]
    ):
        return False
    return True


def get_wt_codings(db: object, trans_flags_dict: dict, trans_id: str) -> tuple:
    """Get matched exons and cds region per transcript id"""

    exon_pos_list = []
    cds_pos_list = []
    cds_frame_list = []
    for feature in db.children(
        trans_id, featuretype=("exon", "CDS"), order_by="start", reverse=False
    ):
        if feature.featuretype == "exon":
            exon_pos_list.append((feature.start, feature.stop))
        elif feature.featuretype == "CDS":
            cds_pos_list.append((feature.start, feature.stop))
            cds_frame_list.append(feature.frame)
    # correct list positions in such a way that the a cds is enclosed
    # by an exon or not all available
    # exons are the scaffold as there can be exons w/o cds but not vv
    cds_pos_list2 = []
    cds_frame_list2 = []
    for exon_start, exon_stop in exon_pos_list:
        exon_has_cds = False
        for cds_counter, (cds_start, cds_stop) in enumerate(cds_pos_list, 0):
            if cds_start >= exon_start and cds_stop <= exon_stop:
                cds_pos_list2.append((cds_start, cds_stop))
                if cds_start > exon_start and cds_stop < exon_stop:
                    trans_flags_dict[trans_id].add(
                        f"cds != exon at no {cds_counter}".format()
                    )
                cds_frame_list2.append(cds_frame_list[cds_counter])
                exon_has_cds = True
        if not exon_has_cds:
            cds_pos_list2.append("NA")
            cds_frame_list2.append("NA")

    return (exon_pos_list, cds_pos_list2, cds_frame_list2)


def merge_features(bp_cds: list, exon_transcripts: list, cds_transcripts: list) -> list:
    """Check if exon is in CDS transcripts for each bp.

    Args:
        exon_transcripts (list): _description_
        cds_transcripts (list): _description_

    Returns:
        list: _description_
    """
    bp_cds_new = []
    for transcript in exon_transcripts:
        if transcript in cds_transcripts:
            bp_cds_new.append(bp_cds[cds_transcripts.index(transcript)])
        else:
            bp_cds_new.append("")
    return bp_cds_new


def get_bp_overlapping_features(db: object, bp_chr: str, bp_pos: int) -> tuple:
    """Get two lists with breakpoint overlapping exons and cds from the database.

    Args:
        db (object): _description_
        bp_chr (str): _description_
        bp_pos (int): _description_

    Returns:
        tuple: _description_
    """

    bp_exons = []
    bp_cds = []
    for efeature in db.region(
        region=f"{bp_chr}:{bp_pos - 1}-{bp_pos + 1}",
        completely_within=False,
        featuretype=("exon", "CDS"),
    ):
        if efeature.featuretype == "exon":
            bp_exons.append(efeature)
        elif efeature.featuretype == "CDS":
            bp_cds.append(efeature)
    # correct positions in such a way that the same position
    # in the exon/CDS list represents the same parental transcript
    # exons are the scaffold as there can be exons w/o cds but not vice versa
    exon_transcripts = [x.attributes["Parent"][0] for x in bp_exons]
    cds_transcripts = [x.attributes["Parent"][0] for x in bp_cds]
    bp_cds = merge_features(bp_cds, exon_transcripts, cds_transcripts)
    return (bp_exons, bp_cds)


def get_parents(db: object, efeature: object) -> tuple:
    """Get transcript and gene parents from exon feature"""
    trans_id = ""
    trans_biotype = ""
    gene_id = ""
    gene_name = ""
    gene_biotype = ""
    description = ""
    count_parents = 0
    for count_parents, parent in enumerate(db.parents(efeature.id), 1):
        if parent.id.startswith("transcript:"):
            trans_id = parent.id
            trans_biotype = parent.attributes["biotype"][0]
        elif parent.id.startswith("gene:"):
            gene_id = parent.attributes.get("gene_id", "")
            gene_name = parent.attributes.get("Name", gene_id)[0]
            gene_biotype = parent.attributes["biotype"][0]
            try:
                description = parent.attributes["description"][0].replace(";", " ")
            except KeyError:
                description = "NA"
    if count_parents > 2:
        pass
        # logger.warning(
        #     "%s has more than two parents. Please check that gene_id %s and \
        #         trans_id %s are correct for the exon at position %s-%s.",
        #         efeature.id, gene_id, trans_id, efeature.start, efeature.stop
        # )
    return (trans_id, trans_biotype, gene_name, gene_biotype, description)
