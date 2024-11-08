"""
This module provides methods to parse the exon/CDS boundaries annotation.
"""

def get_boundary(bp_pos: int, efeature: object) -> str:
    """Check whether the breakpoint position is on an exon or CDS boundary"""
    if not efeature:
        return "NA"
    if bp_pos == efeature.start:
        return "left_boundary"
    if bp_pos == efeature.stop:
        return "right_boundary"
    if efeature.start < bp_pos < efeature.stop:
        return "within"
    return "outside"


def get_boundary2(boundary1: str, boundary2: str, strand1: str, strand2: str) -> str:
    """Check whether boundaries of both fusion partners are in the correct orientation"""

    is_on_boundary = any("left_boundary", "right_boundary")
    stranded = any("+", "-")

    boundary_dict = {
        ("+", "+", "right_boundary", "left_boundary"): "both",
        ("+", "-", "right_boundary", "right_boundary"): "both",
        ("-", "+", "left_boundary", "left_boundary"): "both",
        ("-", "-", "left_boundary", "right_boundary"): "both",
        (stranded, stranded, is_on_boundary, not is_on_boundary): "5prime",
        (stranded, stranded, not is_on_boundary, is_on_boundary): "3prime"
    }
    return boundary_dict.get((strand1, strand2, boundary1, boundary2), "no_match")
