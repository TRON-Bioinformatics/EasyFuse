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

    def __init__(self, cds_id, cds_start, cds_stop, cds_frame, cds_transcript_id):
        self.id = cds_id
        self.start = cds_start
        self.stop = cds_stop
        self.frame = cds_frame
        self.parent_id = cds_transcript_id


    def __repr__(self):
        return f"CDS(id={repr(self.id)}, " \
               f"start={repr(self.start)}, " \
               f"stop={repr(self.stop)}, " \
               f"frame={repr(self.frame)}, " \
               f"parent_id={repr(self.parent_id)})"


    def __eq__(self, other):
        return (
            self.id == other.id and
            self.start == other.start and
            self.stop == other.stop and
            self.frame == other.frame and
            self.parent_id == other.parent_id
        )


    def get_length(self) -> int:
        """Return the CDS length"""
        return self.stop - self.start + 1
