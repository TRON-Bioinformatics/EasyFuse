"""
Module for the CDS class.
"""

class CDS:
    """Class to store CDS information."""
    id: int
    start: int
    stop: int
    frame: int
    parent_id: int

    def __init__(self, cds_id, start, stop, frame, transcript_id):
        self.id = cds_id
        self.start = start
        self.stop = stop
        self.frame = frame
        self.transcript_id = transcript_id


    def __repr__(self):
        return f"CDS(cds_id={repr(self.id)}, " \
               f"start={repr(self.start)}, " \
               f"stop={repr(self.stop)}, " \
               f"frame={repr(self.frame)}, " \
               f"transcript_id={repr(self.transcript_id)})"


    def __eq__(self, other):
        return (
            self.id == other.id and
            self.start == other.start and
            self.stop == other.stop and
            self.frame == other.frame and
            self.transcript_id == other.transcript_id
        )


    def __hash__(self):
        return hash(repr(self))


    def get_length(self) -> int:
        """Return the CDS length"""
        return self.stop - self.start + 1
