"""
write everything to a couple of files
1) context.csv.debug:
all fusions -> used for debugging
2) context.csv:
filtered fusions -> contains only downstream required columns
3) context_seqs.csv.fasta:
contain the fasta sequences of the ft/wt1/wt2 context sequences,
the full length reconstructed transcripts of the gene fusion and
the full length translated peptide sequence of the gene fusion
4) context_seqs.csv_transcripts.fasta (for debugging)
5) context_seqs.csv_peptide.fasta (for debugging)
6) context_seqs.csv.fasta.info :
Length and count of fasta sequences in context_seqs.csv.fasta

"""

import csv

# pylint: disable=E0401
from Bio import SeqIO # type: ignore
from Bio.SeqRecord import SeqRecord # type: ignore

from .file_headers import SHORT_ANNOTATION_HEADER
from .file_headers import FULL_ANNOTATION_HEADER

class OutputHandler:
    """This class handles the output of the fusion annotation tool."""

    def __init__(self, results: list, output_prefix: str, tsl_filter_level: str):
        self.results = results
        self.output_prefix = output_prefix
        self.tsl_filter_level = tsl_filter_level


    def filter_based_on_tsl(self, result: dict) -> bool:
        """Filter fusion events based on TSL level.

        Args:
            result (dict): Dictionary with information on the fusion event

        Returns:
            bool: True if the fusion event should be filtered out
        """
        if (
            result["annotation_bias"]
            or result["wt1_TSL"]
            in self.tsl_filter_level.split(",")
            or result["wt2_TSL"]
            in self.tsl_filter_level.split(",")
        ):
            return True
        return False


    def write_filtered_table(self):
        """Write the filtered annotation table to disk."""
        with open(self.output_prefix, "w", encoding="utf8") as outfile:
            writer = csv.DictWriter(
                outfile,
                fieldnames=SHORT_ANNOTATION_HEADER,
                extrasaction='ignore',
                delimiter=";"
            )
            writer.writeheader()
            for result_dict in self.results:
                if self.filter_based_on_tsl(result_dict):
                    continue
                writer.writerow(result_dict)


    def write_full_table(self):
        """Write the full annotation table to disk."""
        with open(f"{self.output_prefix}.debug", "w", encoding="utf8") as outfile_debug:
            writer = csv.DictWriter(
                outfile_debug,
                fieldnames=FULL_ANNOTATION_HEADER,
                delimiter=";"
            )
            writer.writeheader()
            for result_dict in self.results:
                writer.writerow(result_dict)


    def write_fasta(self):
        """Write the context sequences to a fasta file."""
        with open(f"{self.output_prefix}.fasta", "w", encoding="utf8") as cfasta:
            for result_line in self.results:
                if self.filter_based_on_tsl(result_line):
                    continue
                ftid = result_line["FTID"]
                ctx_id = result_line["context_sequence_id"]
                ctx_bp = result_line["context_sequence_bp"]
                ctx_wt_bp = result_line["context_sequence_wt1_bp"]
                ctx_wt2_bp = result_line["context_sequence_wt2_bp"]
                SeqIO.write(
                    SeqRecord(
                        result_line["context_sequence"],
                        id=f"{ftid}_{ctx_id}_{ctx_bp}_ft",
                        name="",
                        description="",
                    ),
                    cfasta,
                    "fasta",
                )
                SeqIO.write(
                    SeqRecord(
                        result_line["context_sequence_wt1"],
                        id=f"{ftid}_{ctx_id}_{ctx_wt_bp}_wt1",
                        name="",
                        description="",
                    ),
                    cfasta,
                    "fasta",
                )
                SeqIO.write(
                    SeqRecord(
                        result_line["context_sequence_wt2"],
                        id=f"{ftid}_{ctx_id}_{ctx_wt2_bp}_wt2",
                        name="",
                        description="",
                    ),
                    cfasta,
                    "fasta",
                )


    def write_transcript_fasta(self):
        """Write the full length transcripts to a fasta file."""
        with open(f"{self.output_prefix}_transcript.fasta", "w", encoding="utf8") as tfasta:
            for result_line in self.results:
                ftid = result_line["FTID"]
                bp_in_fusion_nt = len(result_line["ft1_cds_transcripts"])
                SeqIO.write(
                    SeqRecord(
                        result_line["fusion_transcript"],
                        id=f"{ftid}_ft_{bp_in_fusion_nt}",
                        name="",
                        description="",
                    ),
                    tfasta,
                    "fasta",
                )


    def write_peptide_fasta(self):
        """Write the full length peptides to a fasta file."""
        with open(f"{self.output_prefix}_peptide.fasta", "w", encoding="utf8") as pfasta:
            for result_line in self.results:
                ftid = result_line["FTID"]
                fusion_type = result_line["type"]
                bp_in_fusion_nt = len(result_line["ft1_cds_transcripts"])
                translation_shift1 = result_line["wt1_frame_at_start"]
                bp_in_fusion_aa = (
                    (bp_in_fusion_nt - translation_shift1) * 10.0 // 3
                ) / 10.0
                SeqIO.write(
                    SeqRecord(
                        result_line["fusion_peptide"],
                        id=f"{ftid}_{bp_in_fusion_aa}_{fusion_type}",
                        name="",
                        description="",
                    ),
                    pfasta,
                    "fasta",
                )


    def write_fasta_info(self):
        """Write the length and count of the fasta sequences to a file."""
        count_context_seqs = 0
        summed_context_seqs_len = 0
        for result_line in self.results:
            if self.filter_based_on_tsl(result_line):
                continue
            count_context_seqs += 3
            summed_context_seqs_len += sum(
                map(
                    len,
                    [
                        result_line["context_sequence"],
                        result_line["context_sequence_wt1"],
                        result_line["context_sequence_wt2"],
                    ],
                )
            )

        with open(f"{self.output_prefix}.fasta.info", "w", encoding="utf8") as info:
            info.write(f"count_context_seqs={count_context_seqs}\n")
            info.write(f"summed_context_seqs_len={summed_context_seqs_len}\n")
