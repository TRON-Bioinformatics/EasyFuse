#!/usr/bin/env python3

"""
@author: TRON (PASO)
@version: 20240129
"""
from argparse import ArgumentParser
import csv
import os
import os.path
import subprocess

from logzero import logger
from count_input_reads import get_input_read_count_from_star


class FusionSummary(object):
    """Collect stats of the run and write them to file"""

    def __init__(
            self,
            input_fusions,
            input_fusion_context_seqs,
            input_requant_counts,
            input_read_stats,
            output_table,
            fusion_tools: list
    ):

        # Check if input files exist
        assert os.path.exists(input_fusions) and os.path.isfile(input_fusions), \
            "Fusion file not found: {}".format(input_fusions)
        self.input_fusions = input_fusions
        assert os.path.exists(input_fusion_context_seqs) and os.path.isfile(input_fusion_context_seqs), \
            "Fusion context seq file not found: {}".format(input_fusion_context_seqs)
        self.input_fusion_context_seqs = input_fusion_context_seqs
        assert os.path.exists(input_requant_counts) and os.path.isfile(input_requant_counts), \
            "Requant counts file not found: {}".format(input_requant_counts)
        self.input_requant_counts = input_requant_counts
        assert os.path.exists(input_read_stats) and os.path.isfile(input_read_stats), \
            "Read stats file not found: {}".format(input_read_stats)
        self.input_read_stats = input_read_stats

        # load aligned read counts
        self.input_read_count = get_input_read_count_from_star(input_read_stats)
        self.output_table = output_table
        self.fusion_tools = fusion_tools
        self.data = {}


    def normalize_counts_cpm(self, count):
        """Returns the CPM of a read count based on the number of original input reads"""
        if self.input_read_count == 0:
            return count
        return str(int(count) / float(self.input_read_count) * 1000000.0)


    def load_context_seqs(self):
        """Parse the context_seq table and load it into a dictionary"""
        out_cols = (
            "BPID",
            "Fusion_Gene",
            "Breakpoint1",
            "Breakpoint2",
            "FTID", # <gene_a>_<bp_a>_<trans_a>_<gene_b>_<bp_b>_<trans_b>
            "context_sequence_id", # md5_hashsum
            "context_sequence_100_id",
            "type",
            "exon_nr",
            "exon_starts",
            "exon_ends",
            "exon_boundary1",
            "exon_boundary2",
            "exon_boundary",
            "bp1_frame",
            "bp2_frame",
            "frame",
            "context_sequence",
            "context_sequence_bp",
            "neo_peptide_sequence",
            "neo_peptide_sequence_bp"
        )
        data = {}
        with open(self.input_fusion_context_seqs) as csvfile:
            csv_reader = csv.DictReader(csvfile, delimiter=';')
            for row in csv_reader:
                # Switch to more unique combination as ID
                ftid = row['FTID']
                if ftid in data:
                    logger.info("FTID already in use")
                data[ftid] = row
        logger.info("Loaded annotation data.")
        return data


    def load_detected_fusions(self):
        """Parse the detected fusions table and load it into a dictionary"""
        out_cols = (
            "BPID", # <chr_a>:<pos_a>:<strand_a>_<chr_b>:<pos_b>:<strand_b>
            "Fusion_Gene", # <gene_a>_<gene_b>
            "Breakpoint1", # <chr_a>:<pos_a>:<strand_a>
            "Breakpoint2", # <chr_b>:<pos_b>:<strand_b>
            "Junction_Reads", # <int>
            "Spanning_Reads", # <int>
            "Sample", # <str>
            "Tool" # <str>
        )
        
        data = {}
        with open(self.input_fusions) as csvfile:
            csv_reader = csv.DictReader(csvfile, delimiter=';')
            for row in csv_reader:
                bpid = row['BPID']
                if bpid not in data:
                    data[bpid] = {}
                data[bpid][row['Tool']] = {
                    'junc': row['Junction_Reads'], 
                    'span': row['Spanning_Reads']
                }
        logger.info("Loaded fusion tool data.")
        return data


    def load_requant_counts(self):
        """Parse the detected fusions table and load it into a dictionary"""
        out_cols = (
            "name", # <FTID>_<context_sequence_id>_<pos>_<ft|wt1|wt2>
            "pos",
            "junc",
            "span",
            "anch",
            "a",
            "b"
        )

        data = {}
        # Columns a and b are not relevant, can be excluded in further merging!
        with open(self.input_requant_counts) as csvfile:
            csv_reader = csv.DictReader(csvfile, delimiter='\t')
            for row in csv_reader:
                ftidplus = row['name']
                id_split = ftidplus.rsplit("_", 3)
                context_sequence_id = id_split[1]
                seq_type = id_split[3]
                if not context_sequence_id in data:
                    data[context_sequence_id] = {}
                if not seq_type in data[context_sequence_id]:
                    data[context_sequence_id][seq_type] = row
        logger.info("Loaded requantification data.")
        return data


    def load_data(self):
        """Join the three data tables context_seq, detected_fusions and requantification"""
        logger.info("Loading input tables")
        context_seqs_data = self.load_context_seqs()
        detected_fusions_data = self.load_detected_fusions()
        requant_data = self.load_requant_counts()
        
        return context_seqs_data, detected_fusions_data, requant_data


    def add_context_data(self, context_seqs_data):
        """Adds the values from fusion annotation to the data"""
        for key in context_seqs_data:
            self.data[key] = {}
            for ele in context_seqs_data[key]:
                self.data[key][ele] = context_seqs_data[key][ele]
        logger.info("Added annotation data.")


    def add_tool_counts(self, detected_fusions_data):
        """Adds the counts from the individual tools to the data"""
        for key in self.data:
            for bpid in detected_fusions_data:
                if self.data[key]["BPID"] == bpid:
                    for tool in self.fusion_tools:
                        if tool in detected_fusions_data[bpid]:
                            self.data[key]["{}_detected".format(tool.lower())] = "1"
                            self.data[key]["{}_junc".format(tool.lower())] = self.normalize_counts_cpm(detected_fusions_data[bpid][tool.lower()]["junc"])
                            self.data[key]["{}_span".format(tool.lower())] = self.normalize_counts_cpm(detected_fusions_data[bpid][tool.lower()]["span"])
                        else:
                            self.data[key]["{}_detected".format(tool.lower())] = "0"
                            self.data[key]["{}_junc".format(tool.lower())] = "0"
                            self.data[key]["{}_span".format(tool.lower())] = "0"
                self.data[key]["tool_count"] = str(len(detected_fusions_data[bpid]))
                self.data[key]["tool_frac"] = str(float(len(detected_fusions_data[bpid]))/len(self.fusion_tools))
        logger.info("Added tool count data.")


    def add_requant_counts(self, requant_data):
        """Adds the counts from requantification to the data"""
        # TODO: Remove redundancy
        for key in self.data:
            ctx_id = self.data[key]['context_sequence_id']
            for seq_type in ("ft", "wt1", "wt2"):
                self.data[key]["{}_bp".format(seq_type)] = requant_data[ctx_id][seq_type]["pos"]
                self.data[key]["{}_junc".format(seq_type)] = self.normalize_counts_cpm(requant_data[ctx_id][seq_type]["junc"])
                self.data[key]["{}_span".format(seq_type)] = self.normalize_counts_cpm(requant_data[ctx_id][seq_type]["span"])
                self.data[key]["{}_anch".format(seq_type)] = requant_data[ctx_id][seq_type]["anch"]
                self.data[key]["{}_bp_cnt".format(seq_type)] = requant_data[ctx_id][seq_type]["pos"]
                self.data[key]["{}_junc_cnt".format(seq_type)] = requant_data[ctx_id][seq_type]["junc"]
                self.data[key]["{}_span_cnt".format(seq_type)] = requant_data[ctx_id][seq_type]["span"]
                self.data[key]["{}_anch_cnt".format(seq_type)] = requant_data[ctx_id][seq_type]["anch"]
        logger.info("Added requantification data.")


    def write_fusion_table(self):
        """Writes the merged data to the final output file"""
        # List columns for final output file
        out_cols = [
            "BPID", # context_seqs + detected_fusions
            "context_sequence_id", # context_seqs
            "FTID", # context_seqs + requant
            "Fusion_Gene", # context_seqs
            "Breakpoint1", # context_seqs + detected_fusions
            "Breakpoint2", # context_seqs + detected_fusions
            "context_sequence_100_id", # context_seqs
            "type", # context_seqs
            "exon_nr", # context_seqs
            "exon_starts", # context_seqs
            "exon_ends", # context_seqs
            "exon_boundary1", # context_seqs
            "exon_boundary2", # context_seqs
            "exon_boundary", # context_seqs
            "bp1_frame", # context_seqs
            "bp2_frame", # context_seqs
            "frame", # context_seqs
            "context_sequence", # context_seqs
            "context_sequence_bp", # context_seqs
            "neo_peptide_sequence", # context_seqs
            "neo_peptide_sequence_bp", # context_seqs
        ]
        # detected_fusions
        for tool in self.fusion_tools:
            out_cols.append("{}_detected".format(tool.lower()))
            out_cols.append("{}_junc".format(tool.lower()))
            out_cols.append("{}_span".format(tool.lower()))
        out_cols.extend([
            "ft_bp", # requant
            "ft_junc", # requant
            "ft_span", # requant
            "ft_anch", # requant
            "wt1_bp", # requant
            "wt1_junc", # requant
            "wt1_span", # requant
            "wt1_anch", # requant
            "wt2_bp", # requant
            "wt2_junc", # requant
            "wt2_span", # requant
            "wt2_anch", # requant
            "ft_bp_cnt", # requant, necessary?
            "ft_junc_cnt", # requant
            "ft_span_cnt", # requant
            "ft_anch_cnt", # requant
            "wt1_bp_cnt", # requant, necessary?
            "wt1_junc_cnt", # requant
            "wt1_span_cnt", # requant
            "wt1_anch_cnt", # requant
            "wt2_bp_cnt", # requant, necessary?
            "wt2_junc_cnt", # requant
            "wt2_span_cnt", # requant
            "wt2_anch_cnt" # requant
        ])

        with open(self.output_table, "w") as outf:
            n = 1
            outf.write(";".join(out_cols) + "\n")
            for key in self.data:
                vals = [self.data[key][ele] if ele in self.data[key] else "NA" for ele in out_cols]
                outf.write(";".join(vals) + "\n")
                n += 1
        return (n, len(out_cols))


    def run(self):
        """Execute individual methods"""

        # load input tables
        context_seqs_data, detected_fusions_data, requant_data = self.load_data()


        self.add_context_data(context_seqs_data)

        self.add_tool_counts(detected_fusions_data)

        self.add_requant_counts(requant_data)

        num_fusions, num_cols = self.write_fusion_table()

        logger.info(
            "Done. Created {} [rows={}, cols={}].".format(
                self.output_table,
                num_fusions,
                num_cols
            )
        )


def main():
    parser = ArgumentParser(description="Merges results to final output table")

    parser.add_argument(
        "--input-fusions",
        dest="input_fusions",
        required=True,
        help="Path to input file with fusions",
    )
    parser.add_argument(
        "--input-fusion-context-seqs",
        dest="input_fusion_context_seqs",
        required=True,
        help="Path to input file with fusion context sequences",
    )
    parser.add_argument(
        "--input-requant-counts",
        dest="input_requant_counts",
        required=True,
        help="Path to input file with requant counts",
    )
    parser.add_argument(
        "--input-read-stats",
        dest="input_read_stats",
        required=True,
        help="Path to input file with read stats",
    )
    parser.add_argument(
        "-o",
        "--output-table",
        dest="output_table",
        help="Specify the merged output file."
    )
    parser.add_argument(
        "--fusion-tools",
        dest="fusion_tools",
        help="Fusion tools."
    )
    args = parser.parse_args()

    summary = FusionSummary(
        input_fusions=args.input_fusions,
        input_fusion_context_seqs=args.input_fusion_context_seqs,
        input_requant_counts=args.input_requant_counts,
        input_read_stats=args.input_read_stats,
        output_table=args.output_table,
        fusion_tools=args.fusion_tools.split(","),
    )
    summary.run()

if __name__ == "__main__":
    main()