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

from .file_headers import FULL_ANNOTATION_HEADER

class OutputHandler:
    """This class handles the output of the fusion annotation tool."""

    def __init__(self, results: list, output_prefix: str, tsl_filter_level: str):
        self.results = results
        self.output_prefix = output_prefix
        self.tsl_filter_level = tsl_filter_level
        self.full_annotation_header_str = ";".join(FULL_ANNOTATION_HEADER)
        self.short_header = FULL_ANNOTATION_HEADER[0:25]
        self.short_header_str = ";".join(self.short_header)
        header_idx = range(len(FULL_ANNOTATION_HEADER))
        self.header_dict = dict(zip(FULL_ANNOTATION_HEADER, header_idx))


    def filter_based_on_tsl(self, result_line) -> bool:
        if (
            result_line["annotation_bias"]
            or result_line["wt1_TSL"]
            in self.tsl_filter_level.split(",")
            or result_line["wt2_TSL"]
            in self.tsl_filter_level.split(",")
        ):
            return True
        return False

    def write_filtered_table(self):
        with open(self.output_prefix, "w", encoding="utf8") as outfile:
            outfile.write(f"{self.short_header_str}\n")
            for result_line in self.results:
                if self.filter_based_on_tsl(result_line):
                    continue
                data_str = ";".join(
                    [
                        str(result_line[col_name])
                        for col_name in self.short_header
                    ]
                )
                outfile.write(f"{data_str}\n")


    def write_full_table(self):
        with open(f"{self.output_prefix}.debug", "w", encoding="utf8") as outfile_debug:
            writer = csv.DictWriter(outfile_debug, fieldnames=FULL_ANNOTATION_HEADER, delimiter=";")
            #outfile_debug.write(f"{self.full_annotation_header_str}\n")
            writer.writeheader()
            for result_dict in self.results:
                writer.writerow(result_dict)
                # for colname in FULL_ANNOTATION_HEADER:
                #     outfile_debug.write(f"{result_dict[colname]};")
                # #res_str = ";".join(map(str, result_line))
                # outfile_debug.write(f"{res_str}\n")


    def write_fasta(self):
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
