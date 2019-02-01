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
import math
import glob
import os
import sys
import pandas as pd
from Bio import pairwise2 # Bio is not available where I run pylint => pylint: disable=E0401
import misc.io_methods as IOMethods
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
        self.chars = ["A", "C", "G", "T"]
        self.overwrite = overwrite
        pd.options.mode.chained_assignment = None
        self.input_read_count = 0

    @staticmethod
    def read_assembly_files(blast_assembly_dir):
        """Read de novo assembled fusion transcripts per gene fusion pair"""
        assembly_data = {}
        assembly_file_list = glob.glob(os.path.join(blast_assembly_dir, "*.blast.out.fusionTranscript"))
        for assembly_file in assembly_file_list:
            fusid = IOMethods.path_leaf(assembly_file).split(".")[0]
            assembly_data[fusid] = {}
            with open(assembly_file, "r") as afile:
                tmp_header = ""
                for i, line in enumerate(afile, 1):
                    if i % 4 == 0:
                        tmp_header = line.split(" ")[0]
                        assembly_data[fusid][tmp_header] = ["NA", "NA", "NA"]
                    elif (i - 1) % 4 == 0:
                        assembly_data[fusid][tmp_header][0] = line.rstrip("\n")
                    elif (i - 2) % 4 == 0:
                        assembly_data[fusid][tmp_header][1] = line.rstrip("\n")
                    else:
                        assembly_data[fusid][tmp_header][2] = line.rstrip("\n")
        return assembly_data

    def entropy_calculation(self, sequence):
        """Return the shannon entropy of a sequence"""
        char_list = list(sequence)
        adj_seq_len = len(sequence)#-char_list.count("N")
        shannon_entropy = 0
        for base in self.chars:
            char_pc = float(char_list.count(base)) / float(adj_seq_len)
            if char_pc != 0:
                shannon_entropy += char_pc * math.log(char_pc, 2)
        return (-1) * shannon_entropy

    @staticmethod
    def similarity_calculation(sequence1, sequence2):
        """Return the score (relative to the max score) of a simple pairwise alignment with 1/-1/-2/-2 scores"""
        alignment = pairwise2.align.globalms(sequence1, sequence2, 1, -1, -2, -2)[0]
        return alignment[2] / float(max(len(sequence1), len(sequence2)))

    @staticmethod
    def filter_annotation_problems(pddf):
        """Identification of annotation problems probably as a result of different annotation versions"""
        count_put_miss_annot = 0
        for i in pddf.index:
            try:
                fusg1, fusg2 = pddf["Fusion_Gene"][i].split("_")
                if not fusg1 in pddf["FTID"][i].upper() or not fusg2 in pddf["FTID"][i].upper():
                    count_put_miss_annot += 1
                    pddf = pddf.drop(i)
            except ValueError:
                count_put_miss_annot += 1
                pddf = pddf.drop(i)
        print("{} putative mis-annotations identified and dropped.".format(count_put_miss_annot))
        return pddf

    def append_to_context_seqs(self, pddf):
        """Append some metrics to the context seq data frame"""
        # urla: the dataframe will require further filtering to be implemented which will be easy starting with a pd df.
        #       furthermore, it garantues to generate a valid csv file.
        #win test:C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\innoData\\fusionProcessing_optV2fltr\\Context_Seqs.csv
        # C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\PRJNA252360_data\\100bp_filtered\\SRR1659951\\classification.tdt
        # C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\PRJNA252360_data\\100bp_filtered\\SRR1659951\\Context_Seqs.csv
        pddf["left_bp_entropy"] = ""
        pddf["right_bp_entropy"] = ""
        pddf["bp_homology_score"] = ""
        pddf["ft_seq_stability"] = 0.0
        for i in pddf.index:
            bp1_sequence = pddf.loc[i, "context_sequence"][pddf.loc[i, "context_sequence_bp"]:]
            bp2_sequence = pddf.loc[i, "context_sequence"][:pddf.loc[i, "context_sequence_bp"]]
            pddf.loc[i, "left_bp_entropy"] = self.entropy_calculation(bp1_sequence)
            pddf.loc[i, "right_bp_entropy"] = self.entropy_calculation(bp2_sequence)
            pddf.loc[i, "bp_homology_score"] = self.similarity_calculation(bp1_sequence, bp2_sequence)
            pddf.loc[i, "ft_seq_stability"] = self.stability_calculation(pddf.loc[i, "context_sequence"], "")
        return pddf

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
                #context_data.loc[i, "{}_detected".format(tool)] = 0
                #context_data.loc[i, "{}_junc".format(tool)] = "NA"
                #context_data.loc[i, "{}_span".format(tool)] = "NA"
                for j in tmp_slice.index:
                    if tool == tmp_slice.loc[j, "Tool"]:
                        context_data.loc[i, "{}_detected".format(tool.lower())] = 1
                        context_data.loc[i, "{}_junc".format(tool.lower())] = self.normalize_counts_cpm(tmp_slice.loc[j, "Junction_Reads"])
                        context_data.loc[i, "{}_span".format(tool.lower())] = self.normalize_counts_cpm(tmp_slice.loc[j, "Spanning_Reads"])
                        context_data.loc[i, "tool_count"] += 1
            context_data.loc[i, "tool_frac"] = float(context_data.loc[i, "tool_count"]) / float(len(fusion_tool_list))
        return context_data

    @staticmethod
    def append_existing_data_hits(context_data, cosmic_data, hbm2_data, icam_high_data, icam_low_data):
        """Add T/F column to indicate whether the respective fusion has been found in other data so far"""
        # so far, the following data is supported:
        # cosmic - https://cancer.sanger.ac.uk/cosmic/fusion?genome=38
        # hbm2 - fusions predicted with easyfuse in the illumina human body map v2 dataset
        # icam high/low - fusions predicted with easyfuse in icam patients with high (in both replicates) or low (in single replicate) confidence
        #                 high confidence is actually a supset of low confidence, hence the name low confidence is actually not correct
        context_data["cosmic_hit"] = False
        context_data["hbm2_hit"] = False
        context_data["icam_high_hit"] = False
        context_data["icam_low_hit"] = False
        for i in context_data.index:
            fg1, fg2 = context_data.loc[i, "Fusion_Gene"].split("_")
            # we check whether fusion gene 1 and 2 are in the respective "names" field and not if the exact same name is matching
            # in order make the search more flexible and easier applicable to different input datasets

            # check cosmic data
            for j in cosmic_data.index:
                if fg1 in cosmic_data.loc[j, "Genes"] and fg2 in cosmic_data.loc[j, "Genes"]:
                    context_data.loc[j, "cosmic_hit"] = True
                    break
            # check hbm2 data
            for j in hbm2_data.index:
                if fg1 in hbm2_data.loc[j, "Genes"] and fg2 in hbm2_data.loc[j, "Genes"]:
                    context_data.loc[j, "hbm2_hit"] = True
                    break
            # check icam high data
            for j in icam_high_data.index:
                if fg1 in icam_high_data.loc[j, "Fusion_Gene"] and fg2 in icam_high_data.loc[j, "Fusion_Gene"]:
                    context_data.loc[j, "icam_high_hit"] = True
                    break
            # check icam low data
            for j in icam_low_data.index:
                if fg1 in icam_low_data.loc[j, "Fusion_Gene"] and fg2 in icam_low_data.loc[j, "Fusion_Gene"]:
                    context_data.loc[j, "icam_low_hit"] = True
                    break
        return context_data

    @staticmethod
    def read_cnt_filtering(pddf):
        """Keep rows fulfilling the defined criteria and drop the rest"""
        for i in pddf.index:
            count_checks = [0, 0, 0, 0, 0, 0, 0, 0]
            # if fusion is supported in sample 1 or 2
            if pddf["junction_ft_1"][i] > 0 and pddf["spanning_ft_1"][i] > 0:
                count_checks[0] = 1
            if pddf["junction_ft_2"][i] > 0 and pddf["spanning_ft_2"][i] > 0:
                count_checks[0] *= 3
                count_checks[1] = 1

            # if wt background is not supported in sample 1 or 2
            if pddf["junction_wt_1"][i] == 0 and pddf["spanning_wt_1"][i] == 0:
                count_checks[2] = 1
            if pddf["junction_wt_2"][i] == 0 and pddf["spanning_wt_2"][i] == 0:
                count_checks[2] *= 3
                count_checks[3] = 1

            # if fusion is support by a long (>= 20) anchor read while wt is not
            if pddf["longest_anchor_ft_1"][i] >= 20 and pddf["longest_anchor_wt_1"][i] < 20:
                count_checks[4] = 1
            if pddf["longest_anchor_ft_2"][i] >= 20 and pddf["longest_anchor_wt_2"][i] < 20:
                count_checks[5] = 1

            # if fusion is on exon boundary
            if pddf["exon_boundary"][i] == "both":
                count_checks[6] = 2
            elif pddf["exon_boundary"][i] == "none":
                count_checks[6] = 0
            else:
                count_checks[6] = -2

            # if fusion transcript is out of frame
            if pddf["frame"][i] == "no_frame":
                count_checks[7] = -99
            else:
                count_checks[7] = 2

            if sum(count_checks) < 10:
                pddf = pddf.drop(i)

        return pddf

    @staticmethod
    def stability_calculation(sequence, dipeptide_score_matrix):
        """Return the stability index as proposed by Guruprasad et al., 1990"""
        dsm = pd.read_table(dipeptide_score_matrix, header=0, index_col=0)
        stab_idx = 0
        for i in range(len(sequence) - 1):
            stab_idx += dsm[sequence[i+1]][sequence[i]]
        return stab_idx * 10 / float(len(sequence))

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

        print("Loading data for sample {} into memory...".format(sample_id))
        if self.check_files(context_seq_file, False):
            context_data = pd.read_csv(context_seq_file, sep=";")
        # append key for later join with requant data
        context_data["ftid_plus"] = context_data["FTID"] + "_" +  context_data["context_sequence_id"]
        redundant_header = list(context_data)
        # urla - todo: filter_annotation_problems searches fusion records where the fusion gene ids from the tool prediction do not match to the
        #              annotation in the context_seqs.csv; I'm not quite sure whether this is a bug in the prediction tools or in our
        #              pipeline
        # urla - note: for now, it is ok NOT to run the following method in order to avoid unintended data removal. The discrepancy between
        #              both annotations are still present, but will not result in errors in downstream processing.
        #context_data = self.filter_annotation_problems(context_data)

        print("Appending normalized fusion counts from {} to the data table.".format(fusion_tools))
        # read the original input read count and append CPM values from individual fusion prediction tools to context data
        with open(os.path.join(os.path.dirname(requant_fltr_file), "Star_org_input_reads.txt"), "r") as rfile:
            self.input_read_count = int(rfile.next())
        if self.check_files(detect_fusion_file, False):
            detect_data = pd.read_csv(detect_fusion_file, sep=";")
        context_data = self.append_tool_cnts_to_context_file(context_data, detect_data, fusion_tools)

        print("Appending normalized requantification counts from filtered and original mapping to the data table.")
        # perform subsequent joins on ftid_plus
        context_data.set_index("ftid_plus", inplace=True)
        # read and append requantification data to context data
        if self.check_files(requant_fltr_file, False):
            requant_fltr_data = pd.read_table(requant_fltr_file, sep=";")
        context_data = context_data.join(requant_fltr_data.set_index("ftid_plus"), how="left")
        if self.check_files(requant_org_file, False):
            requant_org_data = pd.read_table(requant_org_file, sep=";")
        context_data = context_data.join(requant_org_data.set_index("ftid_plus"), lsuffix="_fltr", rsuffix="_org", how="left")

        return context_data.fillna(0), redundant_header


    def run(self, config, icam_run, model_predictions):
        """run the data concatenation"""
        fusion_tools = config.get("general", "fusiontools").split(",")
        # urla - note: with icam_run set to True, two results from technical replicates are expected as input
        print("Looking for required files...")
        # check whether output already exists
        cols_to_aggregate_on = ["FGID", "context_sequence_id", "FTID"]
        
#        joined_table_1 = None
#        joined_table_2 = None
#        joined_table_12 = None

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
                joined_table_2c = joined_table_2.drop(header_list_1, axis=1)
                # perform an inner join of the two replicate data frames
                joined_table_12 = joined_table_1.join(joined_table_2c, lsuffix="_1", rsuffix="_2", how="inner")
                joined_table_12b = joined_table_12.groupby(by=cols_to_aggregate_on, sort=False, as_index=False)
                joined_table_12b.agg(self.custom_data_aggregation).to_csv("{}_fusRank_12.csv".format(self.output), sep=";", index=False)

                # urla - note: the following can currently be considered as "putative developmental leftovers ^^"
    #            assembl_data1 = self.read_assembly_files(assembly_dir1)
    #            assembl_data2 = self.read_assembly_files(assembly_dir2)

    #            print("Appending entropy and homology calculations")
    #            context_data12 = self.append_to_context_seqs(context_data12)
    #            context_data12.to_csv(self.output, sep = "\t", index=True)

    #            print("Performing read support filtering and preparing de novo assembly targets")
    #            context_data12_filtered = self.read_cnt_filtering(context_data12)

    #            with open(os.path.join(assembly_dir1, "fusion_gene_list.txt"), "w") as fgl1, open(os.path.join(assembly_dir2, "fusion_gene_list.txt"), "w") as fgl2:
    #                for fusion_pair in set(context_data12["Fusion_Gene"]):
    #                    fgl1.write("{}\n".format(fusion_pair.replace("_", "--")))
    #                    fgl2.write("{}\n".format(fusion_pair.replace("_", "--")))

        for table in ["1", "2", "12"]:
            summary_file = "{}_fusRank_{}.csv".format(self.output, table)
            if model_predictions and self.check_files(summary_file, True):
                model_exe = config.get("commands", "model_cmd")
                model_path = config.get("otherFiles", "easyfuse_model")
                # append prediction scores based on pre-calculated model
                cmd_model = "{0} --fusion_summary {1} --model_file {2} --output {3}".format(model_exe, summary_file, model_path, "{}.pModelPred.csv".format(summary_file[:-4]))
                Queueing.submit("", cmd_model.split(" "), "", "", "", "", "", "", "none")
        
#        print(joined_table_1)
#        print(joined_table_2)
#        print(joined_table_12)

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
