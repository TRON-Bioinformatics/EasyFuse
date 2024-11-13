"""
This module provides methods to parse the frame annotation.
"""

# pylint: disable=E0401
from breakpoint import Breakpoint

def get_frame(bp: Breakpoint, cds_pos_list: list, cds_frame_list: list) -> tuple:
    """Get/Calculate the frame at the start of the cds and at the breakpoint"""
    # the frame annotation in ensembl and defined as follows:
    # 0: The next full codon (i.e. 3bp codon) starts 0 bases from the current position
    #   => this is the first (0) base of a codon
    # 1: The next full codon (i.e. 3bp codon) starts 1 base from the current position
    #   => this is the third (2) base of a codon
    # 2: The next full codon (i.e. 3bp codon) starts 2 bases from the current position
    #   => this is the second (1) base of a codon
    frame_idx = [0, 2, 1, 0, 2, 1]
    # if the frame at bp is non-determinable, it is set to -1
    frame_at_start = -1
    frame_at_bp = -1
    # zip pos and frame while excluding NAs
    pos_frame_list = [
        (x, int(y)) for x, y in zip(cds_pos_list, cds_frame_list) if not x == "NA"
    ]
    # find breakpoint cds and get the frame at the breakpoint
    for cds_pos, cds_frame in pos_frame_list:
        tmp_frame_idx = frame_idx[
            frame_idx.index(cds_frame) : frame_idx.index(cds_frame) + 3
        ]
        # on pos strand
        if bp.strand == "+":
            if bp.pos == cds_pos[0]:
                frame_at_bp = cds_frame
            elif cds_pos[0] < bp.pos <= cds_pos[1]:
                frame_at_bp = tmp_frame_idx[(bp.pos - (cds_pos[0] - 1)) % 3]
        # on neg strand
        else:
            if bp.pos == cds_pos[1]:
                frame_at_bp = cds_frame
            elif cds_pos[0] <= bp.pos < cds_pos[1]:
                frame_at_bp = tmp_frame_idx[(cds_pos[1] - (bp.pos - 1)) % 3]
    # get the starting frame of the cds
    if not len(pos_frame_list) == 0:
        if bp.strand == "+":
            frame_at_start = pos_frame_list[0][1]
        else:
            frame_at_start = pos_frame_list[-1][1]
    # return frame at the beginning of the first cds
    # and at the breakpoint position, both corrected according the the reading strand
    return frame_at_start, frame_at_bp
