#!/usr/bin/env python

"""
Combining of all relevant infos into a final table
Calculations of additional metrics
Filtering of the final data table

@author: BNT (URLA)
@version: 20190118
"""

from __future__ import print_function, division
import argparse
import os
import sys
import pandas as pd
import misc.queue as Queueing

# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class DataJoining(object):
    """Select alignments belonging to putative fusions from an s/bam file"""
    def __init__(self, input_dir, id1, id2, output, overwrite):
        """Parameter initialization"""
        self.input_dir = input_dir
        self.id1 = id1
        self.id2 = id2
        self.output = output
        self.overwrite = overwrite
        pd.options.mode.chained_assignment = None
        self.input_read_count = 0
        self.easyfuse_path = os.path.dirname(os.path.realpath(__file__))

    def append_tool_cnts_to_context_file(self, context_data, detected_fusions, fusion_tool_list):
        """Add data to context seq table"""
        # Init new columns
        for tool in fusion_tool_list:
            context_data["{}_detected".format(tool.lower())] = 0
            context_data["{}_junc".format(tool.lower())] = 0
            context_data["{}_span".format(tool.lower())] = 0
        context_data["tool_count"] = 0
        context_data["tool_frac"] = 0.0

        # iterate over context seq data
        for i in context_data.index:
            # get tool prediction for each FGID
            tmp_slice = detected_fusions.loc[detected_fusions["FGID"] == context_data["FGID"][i]]
            # check for each tool, whether a prediction was made by it and update the respective columns
            for tool in fusion_tool_list:
                # pre-set values to 0 or NA indicating no prediction with this tool
                # line wise initialization shouldn't be required as it was initialized with this values in the beginning
                for j in tmp_slice.index:
                    if tool == tmp_slice.loc[j, "Tool"]:
                        context_data.loc[i, "{}_detected".format(tool.lower())] = 1
                        context_data.loc[i, "{}_junc".format(tool.lower())] = self.normalize_counts_cpm(tmp_slice.loc[j, "Junction_Reads"])
                        context_data.loc[i, "{}_span".format(tool.lower())] = self.normalize_counts_cpm(tmp_slice.loc[j, "Spanning_Reads"])
                        context_data.loc[i, "tool_count"] += 1
            context_data.loc[i, "tool_frac"] = float(context_data.loc[i, "tool_count"]) / float(len(fusion_tool_list))
        return context_data

    def normalize_counts_cpm(self, count):
        """Returns the CPM of a read count based on the number of original input reads"""
        if self.input_read_count == 0:
            return count
        return count / self.input_read_count * 1000000

    def create_joined_table(self, sample_id, fusion_tools):
        """Join the three data tables context_seq, detected_fusions and requantification"""
        # define path' to the context seq, detected fusion and re-quantification files
        context_seq_file = os.path.join(self.input_dir, "Sample_{}".format(sample_id), "fetchdata", "fd_1_tool", "fetched_contextseqs", "Context_Seqs.csv")
        detect_fusion_file = os.path.join(self.input_dir, "Sample_{}".format(sample_id), "fetchdata", "fd_1_tool", "fetched_fusions", "Detected_Fusions.csv")
        requant_fltr_file = os.path.join(self.input_dir, "Sample_{}".format(sample_id), "fetchdata", "fd_1_tool", "classification", "classification_fltr.tdt")
        requant_org_file = os.path.join(self.input_dir, "Sample_{}".format(sample_id), "fetchdata", "fd_1_tool", "classification", "classification_org.tdt")
        requant_cnt_fltr_file = os.path.join(self.input_dir, "Sample_{}".format(sample_id), "fetchdata", "fd_1_tool", "classification", "classification_fltr.tdt.counts")
        requant_cnt_org_file = os.path.join(self.input_dir, "Sample_{}".format(sample_id), "fetchdata", "fd_1_tool", "classification", "classification_org.tdt.counts")

        print("Loading data for sample {} into memory...".format(sample_id))
        if self.check_files(context_seq_file, False):
            context_data = pd.read_csv(context_seq_file, sep=";")
        # append key for later join with requant data
        context_data["ftid_plus"] = context_data["FTID"] + "_" +  context_data["context_sequence_id"]
        redundant_header = list(context_data)

        print("Appending normalized fusion counts from {} to the data table.".format(fusion_tools))
        # read the original input read count and append CPM values from individual fusion prediction tools to context data
        with open(os.path.join(os.path.dirname(requant_fltr_file), "Star_org_input_reads.txt"), "r") as rfile:
            self.input_read_count = int(rfile.next())
        if self.check_files(detect_fusion_file, False):
            detect_data = pd.read_csv(detect_fusion_file, sep=";")
        context_data = self.append_tool_cnts_to_context_file(context_data, detect_data, fusion_tools)

        print("Appending (normalized) requantification counts from filtered and original mapping to the data table.")
        # perform subsequent joins on ftid_plus
        context_data.set_index("ftid_plus", inplace=True)
        # read and append requantification data (cpm) to context data
        if self.check_files(requant_fltr_file, False):
            requant_fltr_data = pd.read_csv(requant_fltr_file, sep=";")
        context_data = context_data.join(requant_fltr_data.set_index("ftid_plus"), how="left")
        if self.check_files(requant_org_file, False):
            requant_org_data = pd.read_csv(requant_org_file, sep=";")
        context_data = context_data.join(requant_org_data.set_index("ftid_plus"), lsuffix="_fltr", rsuffix="_org", how="left")
        # read and append requantification data (count) to context data
        if self.check_files(requant_cnt_fltr_file, False):
            requant_cnt_fltr_data = pd.read_csv(requant_cnt_fltr_file, sep=";")
        context_data = context_data.join(requant_cnt_fltr_data.set_index("ftid_plus"), how="left")
        if self.check_files(requant_cnt_org_file, False):
            requant_cnt_org_data = pd.read_csv(requant_cnt_org_file, sep=";")
        context_data = context_data.join(requant_cnt_org_data.set_index("ftid_plus"), lsuffix="_cnt_fltr", rsuffix="_cnt_org", how="left")

        return context_data.fillna(0), redundant_header


    def run(self, config, icam_run, model_predictions):
        """run the data concatenation"""
        fusion_tools = config.get("general", "fusiontools").split(",")
        # urla - note: with icam_run set to True, two results from technical replicates are expected as input
        print("Looking for required files...")
        # check whether output already exists
        cols_to_aggregate_on = ["FGID", "context_sequence_id", "FTID"]
        
        if not self.overwrite and self.check_files("{}_fusRank_1.csv".format(self.output), True):
            print("Found pre-calculated output files. If you would like to re-calculate everything, re-start with the --overwrite parameter supplied.")
            joined_table_1 = pd.read_csv("{}_fusRank_1.csv".format(self.output), sep=";")
            if icam_run and self.check_files("{}_fusRank_2.csv".format(self.output), True) and self.check_files("{}_fusRank_12.csv".format(self.output), True):
                joined_table_2 = pd.read_csv("{}_fusRank_2.csv".format(self.output), sep=";")
                joined_table_12 = pd.read_csv("{}_fusRank_12.csv".format(self.output), sep=";")
        else:
            joined_table_1, header_list_1 = self.create_joined_table(self.id1, fusion_tools)

#            print("Headerlist: {}".format(header_list_1))

            joined_table_1b = joined_table_1.groupby(by=cols_to_aggregate_on, sort=False, as_index=False)
            joined_table_1b.agg(self.custom_data_aggregation).to_csv("{}_fusRank_1.csv".format(self.output), sep=";", index=False)
            if icam_run:
                joined_table_2, _ = self.create_joined_table(self.id2, fusion_tools)
                joined_table_2b = joined_table_2.groupby(by=cols_to_aggregate_on, sort=False, as_index=False)
                joined_table_2b.agg(self.custom_data_aggregation).to_csv("{}_fusRank_2.csv".format(self.output), sep=";", index=False)

                print("Merging data...")
#                joined_table_1 = joined_table_1.drop(cols_to_drop, axis=1)
#                cols_to_drop.extend([i for i in cols_from_one if i not in cols_to_drop])
                # drop columns from context seq in tab 2 as they are identical to the ones from tab 1
                joined_table_2c = joined_table_2.drop(header_list_1, axis=1, errors="ignore")
                # perform an inner join of the two replicate data frames
                joined_table_12 = joined_table_1.join(joined_table_2c, lsuffix="_1", rsuffix="_2", how="inner")
                joined_table_12b = joined_table_12.groupby(by=cols_to_aggregate_on, sort=False, as_index=False)
                joined_table_12b.agg(self.custom_data_aggregation).to_csv("{}_fusRank_12.csv".format(self.output), sep=";", index=False)

        for table in ["1", "2"]:
            summary_file = "{}_fusRank_{}.csv".format(self.output, table)
            if model_predictions and self.check_files(summary_file, True):
                model_path = config.get("otherFiles", "easyfuse_model")
                model_threshold = config.get("general", "model_pred_threshold")
                # append prediction scores based on pre-calculated model
                cmd_model = "{0} --fusion_summary {1} --model_file {2} --prediction_threshold {3} --output {4}".format(os.path.join(self.easyfuse_path, "R", "R_model_prediction.R"), summary_file, model_path, model_threshold, "{}.pModelPred.csv".format(summary_file[:-4]))
                Queueing.submit("", cmd_model.split(" "), "", "", "", "", "", "", "", "", "none")

        if icam_run:
            return (self.count_records(joined_table_1), self.count_records(joined_table_2), self.count_records(joined_table_12))
        return self.count_records(joined_table_1)

    @staticmethod
    def custom_data_aggregation(pd_series):
        """Removes duplicates from a pandas series fields and creates a joined string from the resulting list"""
        # urla - note: this is an adaptation from https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
        #              in a way that it does not return a list but a string to print from the list
        seen = set()
        seen_add = seen.add
        return ";".join(str(e) for e in [x for x in pd_series if not (x in seen or seen_add(x))])

    @staticmethod
    def check_files(file_path, dont_quit):
        """Checks whether a file or dir exists. Returns true or sys-exits with error msg"""
        if not os.path.isfile(file_path) and not os.path.isdir(file_path):
            if dont_quit:
                return False
            else:
                print("Error 99: File or directory not found or just not at the expected location. Looked at {}".format(file_path))
                sys.exit(99)
        else:
            print("Found {}".format(file_path))
        return True

    @staticmethod
    def count_records(pd_data_frame):
        """count different fields"""
        gene_blacklist = ["HLA", "IG", "RP", "LINC"]
        counter = []
#        print("pdx: {}".format(list(pd_data_frame)))
        # how many fusion genes have been predicted
        counter.append(len(set(pd_data_frame["Fusion_Gene"].tolist())))
        # how many fusion genes generate NO no-frame fusion (in/out/neo frame)
        counter.append(len(set(pd_data_frame.loc[pd_data_frame["frame"] != "no_frame"]["Fusion_Gene"].tolist())))
        # how many fusion genes are fused on exon boundaries (both fusion genes must be on the boundary)
        counter.append(len(set(pd_data_frame.loc[pd_data_frame["exon_boundary"] == "both"]["Fusion_Gene"].tolist())))
        # how many fusion genes have a predicted neo peptide sequence
        counter.append(len(set(pd_data_frame.loc[pd_data_frame["neo_peptide_sequence"].notna()]["Fusion_Gene"].tolist())))
        # how many fusion genes do NOT belong to either HLA/IGx/MT/Ribosomal gene families
        counter.append(len(set(pd_data_frame.loc[~pd_data_frame["Fusion_Gene"].str.contains("|".join(gene_blacklist), na=False)]["Fusion_Gene"].tolist())))
        # how many fusion genes have been predicted by at least 2 tools
        #counter.append(len(set(pd_data_frame.loc[pd_data_frame["ToolCount"] >= 2]["Fusion_Gene"].tolist())))
        # how many fusion genes fulfill all of the previously mentioned criteria
        hq_fus_candidates = set(pd_data_frame.loc[(pd_data_frame["frame"] != "no_frame") & (pd_data_frame["exon_boundary"] == "both") & (pd_data_frame["neo_peptide_sequence"].notna()) & (~pd_data_frame["Fusion_Gene"].str.contains("|".join(gene_blacklist), na=False))]["Fusion_Gene"].tolist())
        counter.append(len(hq_fus_candidates))
        return (counter, hq_fus_candidates)

def main():
    """Parse command line arguments and start script"""
    parser = argparse.ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input_dir', dest='input_dir', help='Input root directory containing the sample file directories', required=True)
    parser.add_argument('--id1', dest='id1', help='Identifier of the first sample replicate', required=True)
    parser.add_argument('--id2', dest='id2', help='Identifier of the second sample replicate', required=True)
    parser.add_argument('-o', '--output_prefix', dest='output', help='Specify output file prefix', required=True)
    parser.add_argument('-c', '--config', dest='output', help='Specify output file prefix', required=True)
    parser.add_argument('-f', '--fus_tool', dest='fustool', help='List of fusion tools to consider')
    parser.add_argument('--overwrite', dest='overwrite', help=argparse.SUPPRESS, default=False, action='store_true')
    parser.add_argument('--icam_run', dest='overwrite', help=argparse.SUPPRESS, default=False, action='store_true')
    args = parser.parse_args()

    fusrank = DataJoining(args.input_dir, args.id1, args.id2, args.output, args.overwrite)
    fusrank.run(args.fustool, args.icam_run, args.model_predictions)

if __name__ == '__main__':
    main()
