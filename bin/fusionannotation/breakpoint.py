"""
Module for the Breakpoint class.
"""

# pylint: disable=E0401
from exon_validation import merge_features

class Breakpoint:
    """Class to store breakpoint information"""
    def __init__(self, chrom: str, pos: int, strand: str):
        self.chrom = chrom
        self.pos = pos
        self.strand = strand

    def __str__(self):
        return f"{self.chrom}:{self.pos}:{self.strand}"


    def get_overlapping_features(self, db: object) -> tuple:
        """Get two lists with breakpoint overlapping exons and cds from the database.

        Args:
            db (object): GFF database object

        Returns:
            tuple: Exons and CDS overlapping the breakpoint
        """

        bp_exons = []
        bp_cds = []
        for efeature in db.region(
            region=f"{self.chrom}:{self.pos - 1}-{self.pos + 1}",
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


    def get_boundary(self, efeature: object) -> str:
        """Check whether the breakpoint position is on an exon or CDS boundary"""
        if not efeature:
            return "NA"
        if self.pos == efeature.start:
            return "left_boundary"
        if self.pos == efeature.stop:
            return "right_boundary"
        if efeature.start < self.pos < efeature.stop:
            return "within"
        return "outside"
