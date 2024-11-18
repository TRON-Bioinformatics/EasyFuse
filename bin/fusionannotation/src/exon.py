"""
Module for the Exon class.
"""

class Exon:
    """Class to store exon information."""
    def __init__(self, exon_id, exon_start, exon_stop, exon_parent_id):
        self.id = exon_id
        self.start = exon_start
        self.stop = exon_stop
        self.parent_id = exon_parent_id


    def __repr__(self):
        return f"Exon(id={repr(self.id)}, " \
               f"start={repr(self.start)}, " \
               f"stop={repr(self.stop)}, " \
               f"parent_id={repr(self.parent_id)})"


    def __eq__(self, other):
        return (
            self.id == other.id and
            self.start == other.start and
            self.stop == other.stop and
            self.parent_id == other.parent_id
        )


    def get_length(self) -> int:
        """Return the exon length"""
        return self.stop - self.start + 1
