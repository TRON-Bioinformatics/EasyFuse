"""
Module for storing transcript information.
"""

class Transcript:
    """Class to store transcript information"""

    def __init__(self, transcript_id, transcript_biotype, gene_name, gene_biotype, description):
        self.transcript_id = transcript_id
        self.transcript_biotype = transcript_biotype
        self.gene_name = gene_name
        self.gene_biotype = gene_biotype
        self.description = description
        self.exon_pos_list = []
        self.cds_pos_list = []
        self.frame = None
        self.tsl = None


    def set_exon_pos_list(self, exon_pos_list: list):
        """Sets the exon position list.

        Args:
            exon_pos_list (list): Exon positions
        """
        self.exon_pos_list = exon_pos_list


    def set_cds_pos_list(self, cds_pos_list: list):
        """Sets the cds position list.

        Args:
            cds_pos_list (list): CDS positions
        """
        self.cds_pos_list = cds_pos_list


    def set_frame(self, frame: int):
        """Sets the frame.

        Args:
            frame (int): Frame
        """
        self.frame = frame


    def set_tsl(self, tsl: str):
        """Sets the Transcript Support Level (TSL).

        Args:
            tsl (str): TSL
        """
        self.tsl = tsl


    def __eq__(self, other):
        return (
            self.transcript_id == other.transcript_id and
            self.transcript_biotype == other.transcript_biotype and
            self.gene_name == other.gene_name and
            self.gene_biotype == other.gene_biotype and
            self.description == other.description
        )


    def __repr__(self):
        return f"Transcript(transcript_id={repr(self.transcript_id)}, " \
               f"transcript_biotype={repr(self.transcript_biotype)}, " \
               f"gene_name={repr(self.gene_name)}, " \
               f"gene_biotype={repr(self.gene_biotype)}, " \
               f"description={repr(self.description)}"
