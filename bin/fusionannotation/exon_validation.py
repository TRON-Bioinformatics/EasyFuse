"""
This module performs validation of exons.
"""

# pylint: disable=E0401
from transcript import Transcript

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


def extract_exons_from_transcript(db: object, trans_id: str) -> list:
    """Extract exon features from a transcript id"""
    exon_pos_list = []
    for feature in db.children(
        trans_id, featuretype=("exon"), order_by="start", reverse=False
    ):
        exon_pos_list.append((feature.start, feature.stop))
    return exon_pos_list


def extract_cds_frame_from_transcript(db: object, trans_id: str) -> list:
    """Extract cds frame from a transcript id"""
    cds_pos_list = []
    cds_frame_list = []
    for feature in db.children(
        trans_id, featuretype=("CDS"), order_by="start", reverse=False
    ):
        cds_pos_list.append((feature.start, feature.stop))
        cds_frame_list.append(feature.frame)
    return cds_pos_list, cds_frame_list


def extract_features_from_transcript(db: object, trans_id: str) -> tuple:
    """Extract exon and cds features from a transcript id"""
    exon_pos_list = extract_exons_from_transcript(db, trans_id)
    cds_pos_list, cds_frame_list = extract_cds_frame_from_transcript(db, trans_id)
    return (exon_pos_list, cds_pos_list, cds_frame_list)


def filter_cds(exon_pos: list, cds_pos: list, cds_frame: list, trans_flags_dict: dict, trans_id: str) -> tuple:
    """Get matched exons and cds region per transcript id
    Check if cds is enclosed by exon and if cds != exon"""

    # correct list positions in such a way that the a cds is enclosed
    # by an exon or not all available
    # exons are the scaffold as there can be exons w/o cds but not vv
    cds_pos_filtered = []
    cds_frame_filtered = []
    for exon_start, exon_stop in exon_pos:
        exon_has_cds = False
        for cds_counter, (cds_start, cds_stop) in enumerate(cds_pos, 0):
            if cds_start >= exon_start and cds_stop <= exon_stop:
                cds_pos_filtered.append((cds_start, cds_stop))
                if cds_start > exon_start and cds_stop < exon_stop:
                    trans_flags_dict[trans_id].add(
                        f"cds != exon at no {cds_counter}"
                    )
                cds_frame_filtered.append(cds_frame[cds_counter])
                exon_has_cds = True
        if not exon_has_cds:
            cds_pos_filtered.append("NA")
            cds_frame_filtered.append("NA")

    return (cds_pos_filtered, cds_frame_filtered)


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


def get_parents(db: object, efeature: object) -> Transcript:
    """Get transcript and gene parents from exon feature"""
    trans_id = ""
    trans_biotype = ""
    gene_id = ""
    gene_name = ""
    gene_biotype = ""
    description = ""
    #count_parents = 0
    for parent in db.parents(efeature.id):
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
    #if count_parents > 2:
    #    pass
        # logger.warning(
        #     "%s has more than two parents. Please check that gene_id %s and \
        #         trans_id %s are correct for the exon at position %s-%s.",
        #         efeature.id, gene_id, trans_id, efeature.start, efeature.stop
        # )
    return Transcript(
        trans_id, 
        trans_biotype, 
        gene_name, 
        gene_biotype, 
        description
    )
