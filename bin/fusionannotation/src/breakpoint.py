"""
Module for the Breakpoint class.
"""


class Breakpoint:
    """Class to store breakpoint information"""
    def __init__(self, chrom: str, pos: int, strand: str):
        self.chrom = chrom
        self.pos = pos
        self.strand = strand

    def __repr__(self):
        return f"{self.chrom}:{self.pos}:{self.strand}"


    def __eq__(self, other: object) -> bool:
        return (
            self.chrom == other.chrom
            and self.pos == other.pos
            and self.strand == other.strand
        )


    def __hash__(self):
        return hash(repr(self))


    def get_boundary(self, feature: object) -> str:
        """Check whether the breakpoint position is on an exon or CDS boundary"""
        if not feature:
            return "NA"
        if self.pos == feature.start:
            return "left_boundary"
        if self.pos == feature.stop:
            return "right_boundary"
        if feature.start < self.pos < feature.stop:
            return "within"
        return "outside"


    def get_frame(self, cds_pos_list: list) -> tuple:
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

        pos_frame_list = [
            (cds, cds.frame) for cds in cds_pos_list if cds is not None
        ]
        # find breakpoint cds and get the frame at the breakpoint
        for cds, cds_frame in pos_frame_list:
            tmp_frame_idx = frame_idx[
                frame_idx.index(cds_frame) : frame_idx.index(cds_frame) + 3
            ]
            # on pos strand
            if self.strand == "+":
                if self.pos == cds.start:
                    frame_at_bp = cds_frame
                elif cds.start < self.pos <= cds.stop:
                    frame_at_bp = tmp_frame_idx[(self.pos - (cds.start - 1)) % 3]
            # on neg strand
            else:
                if self.pos == cds.stop:
                    frame_at_bp = cds_frame
                elif cds.start <= self.pos < cds.stop:
                    frame_at_bp = tmp_frame_idx[(cds.stop - (self.pos - 1)) % 3]
        # get the starting frame of the cds
        if not len(pos_frame_list) == 0:
            if self.strand == "+":
                frame_at_start = pos_frame_list[0][1]
            else:
                frame_at_start = pos_frame_list[-1][1]
        # return frame at the beginning of the first cds
        # and at the breakpoint position, both corrected according the the reading strand
        return frame_at_start, frame_at_bp
