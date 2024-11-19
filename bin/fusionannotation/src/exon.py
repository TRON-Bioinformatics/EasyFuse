"""
Module for the Exon class.
"""

class Exon:
    """Class to store exon information."""
    def __init__(self, exon_id, start, stop, transcript_id):
        self.id = exon_id
        self.start = start
        self.stop = stop
        self.transcript_id = transcript_id


    def __repr__(self):
        return f"Exon(exon_id={repr(self.id)}, " \
               f"start={repr(self.start)}, " \
               f"stop={repr(self.stop)}, " \
               f"transcript_id={repr(self.transcript_id)})"


    def __eq__(self, other):
        return (
            self.id == other.id and
            self.start == other.start and
            self.stop == other.stop and
            self.transcript_id == other.transcript_id
        )


    def get_length(self) -> int:
        """Return the exon length"""
        return self.stop - self.start + 1
