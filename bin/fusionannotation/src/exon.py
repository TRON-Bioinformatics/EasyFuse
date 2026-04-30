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


    def __hash__(self):
        return hash(repr(self))


    def get_length(self) -> int:
        """Return the exon length"""
        return self.stop - self.start + 1


def join_exon_starts(exons: list) -> str:
    """_summary_

    Args:
        exons (list): _description_

    Returns:
        str: _description_
    """
    return "|".join([str(exon.start) for exon in exons])


def join_exon_stops(exons: list) -> str:
    """_summary_

    Args:
        exons (list): _description_

    Returns:
        str: _description_
    """
    return "|".join([str(exon.stop) for exon in exons])


def get_exon_ranges(exons: list, neo_pep_until_bp_nuc: int) -> list:
    """_summary_

    Args:
        exons (list): _description_
        neo_pep_until_bp_nuc (int): _description_

    Returns:
        list: _description_
    """
    exon_pos = 1
    summed_len = 0
    for exon_pos, exon_length in enumerate([exon.stop - exon.start for exon in exons], 1):
        summed_len += exon_length
        if neo_pep_until_bp_nuc <= summed_len:
            break
    lookup_exons = exons[:exon_pos]
    exon_starts = join_exon_starts(lookup_exons)
    exon_stops = join_exon_stops(lookup_exons)
    return (exon_starts, exon_stops)


def get_exon_ranges_reverse(exons: list, neo_pep_until_bp_nuc: int) -> tuple:
    """_summary_

    Args:
        exons (list): _description_
        neo_pep_until_bp_nuc (int): _description_

    Returns:
        tuple: _description_
    """
    exon_pos = 1
    summed_len = 0
    for exon_pos, exon_length in enumerate([exon.stop - exon.start for exon in exons][::-1], 1):
        summed_len += exon_length
        if neo_pep_until_bp_nuc <= summed_len:
            break
    lookup_exons = exons[-exon_pos:]
    exon_starts = join_exon_starts(lookup_exons)
    exon_stops = join_exon_stops(lookup_exons)
    return (exon_starts, exon_stops)
