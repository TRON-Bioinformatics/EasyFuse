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
        
        # in comparison to the previous implementation, pre-existing results will always be overwritten as it turned out that the
        # previous implementation is to error prone...
        joined_table_1, header_list_1 = self.create_joined_table(self.id1, fusion_tools)
        joined_table_1b = joined_table_1.groupby(by=cols_to_aggregate_on, sort=False, as_index=False)
        joined_table_1b.agg(self.custom_data_aggregation).to_csv("{}_fusRank_1.csv".format(self.output), sep=";", index=False)
        if icam_run:
            joined_table_2, _ = self.create_joined_table(self.id2, fusion_tools)
            joined_table_2b = joined_table_2.groupby(by=cols_to_aggregate_on, sort=False, as_index=False)
            joined_table_2b.agg(self.custom_data_aggregation).to_csv("{}_fusRank_2.csv".format(self.output), sep=";", index=False)
            print("Merging data...")
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
        # re-read the table with prediction for filter counting
        model_table_1 = pd.read_csv("{}_fusRank_1.pModelPred.csv".format(self.output), sep=";")

        if icam_run:
            model_table_2 = pd.read_csv("{}_fusRank_2.pModelPred.csv".format(self.output), sep=";")
            #model_table_2.drop(header_list_1, axis=1, errors="ignore", inplace=True)
            model_table_12 = model_table_1.join(model_table_2.drop(header_list_1, axis=1, errors="ignore"), lsuffix="_1", rsuffix="_2", how="inner")
            model_table_12b = model_table_12.groupby(by=cols_to_aggregate_on, sort=False, as_index=False)
            model_table_12b.agg(self.custom_data_aggregation).to_csv("{}_fusRank_12.pModelPred.csv".format(self.output), sep=";", index=False)
            model_table_12 = pd.read_csv("{}_fusRank_12.pModelPred.csv".format(self.output), sep=";")
            return (self.count_records(model_table_1, False, "F1"), self.count_records(model_table_2, False, "F2"), self.count_records(model_table_12, True, "F12"))
        return self.count_records(model_table_1, False, "F1N")

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
    def count_records_old(pd_data_frame):
        """ the old count records method """
        #gene_blacklist = ["HLA", "IG", "RP", "LINC"]
        gene_blacklist = ["HLA", "IG"]
        counter = []
#        print("pdx: {}".format(list(pd_data_frame)))
        # how many fusion genes have been predicted
#        counter.append(len(set(pd_data_frame["Fusion_Gene"].tolist())))
        counter.append(pd_data_frame["Fusion_Gene"].nunique())
        # how many fusion genes generate NO no-frame fusion (in/out/neo frame)
#        counter.append(len(set(pd_data_frame.loc[pd_data_frame["frame"] != "no_frame"]["Fusion_Gene"].tolist())))
        counter.append(pd_data_frame[pd_data_frame["frame"] != "no_frame"]["Fusion_Gene"].nunique())
        # how many fusion genes are fused on exon boundaries (both fusion genes must be on the boundary)
        counter.append(pd_data_frame[pd_data_frame["exon_boundary"] == "both"]["Fusion_Gene"].nunique())
        # how many fusion genes have a predicted neo peptide sequence
        counter.append(pd_data_frame[(pd_data_frame["neo_peptide_sequence"].notna()) & (pd_data_frame["neo_peptide_sequence"] != "0")]["Fusion_Gene"].nunique())
        # how many fusion genes do NOT belong to either HLA/IGx/MT/Ribosomal gene families
        counter.append(pd_data_frame[~pd_data_frame["Fusion_Gene"].str.contains("|".join(gene_blacklist), na=False)]["Fusion_Gene"].nunique())
        # how many fusions got at least 1 spanning or junction read on the ft during requanification
        testdf = pd_data_frame[["Fusion_Gene", "ft_junc_cnt_org", "ft_span_cnt_org"]].groupby(by="Fusion_Gene")
        counter.append(((testdf["ft_junc_cnt_org"].max() > 0) | (testdf["ft_span_cnt_org"].max() > 0)).sum())
        # how many fusions with at least one positive transcript according to the random forrest model
        counter.append(int((pd_data_frame[["Fusion_Gene","prediction_prob"]].groupby(['Fusion_Gene']).max() >= 0.75).sum()))
        
        counter.append(pd_data_frame[~pd_data_frame["Fusion_Gene"].str.contains("|".join(gene_blacklist), na=False)]["Fusion_Gene"].nunique())
        counter.append(pd_data_frame[~pd_data_frame["Fusion_Gene"].str.contains("|".join(gene_blacklist), na=False)]["Fusion_Gene"].nunique())
        counter.append(pd_data_frame[~pd_data_frame["Fusion_Gene"].str.contains("|".join(gene_blacklist), na=False)]["Fusion_Gene"].nunique())
        # how many fusion genes have been predicted by at least 2 tools
        #counter.append(len(set(pd_data_frame.loc[pd_data_frame["ToolCount"] >= 2]["Fusion_Gene"].tolist())))
        # how many fusion genes fulfill all of the previously mentioned criteria
        hq_fus_candidates = set(pd_data_frame.loc[(pd_data_frame["frame"] != "no_frame") & (pd_data_frame["exon_boundary"] == "both") & (pd_data_frame["neo_peptide_sequence"].notna()) & (~pd_data_frame["Fusion_Gene"].str.contains("|".join(gene_blacklist), na=False))]["Fusion_Gene"].tolist())
        counter.append(len(hq_fus_candidates))
        return (counter, hq_fus_candidates)
    
    @staticmethod
    def count_records(pd_data_frame, is_merged, name):
        """ stub """
        print("Processing: {}".format(name))
        print("Header: {}".format(list(pd_data_frame)))
        gene_blacklist = ["HLA", "IG"]
        # list of sets to store (un)filtered fusion gene names
        # 0:type; 1:boundary; 2:frame; 3:pepseq; 4:counts; 5:blacklist; 6:prediction;
        # 7:exonError; 8:unfiltered; 9:allFilter; 10:allButPredFilter
        fusion_gene_set_list = []
        df_rep = pd_data_frame
        # try to correct annotation issues
        ftid_dict = {}
        for i in df_rep.index:
            ftid = df_rep.loc[i, "FTID"]
            if ftid not in ftid_dict:
                ftid_dict[ftid] = [df_rep.loc[i, "Fusion_Gene"]]
            else:
                ftid_dict[ftid].append(df_rep.loc[i, "Fusion_Gene"])
        
        for ftid in ftid_dict:
            # check if multiple fusion_gene strings are associated with the ftid
            if len(ftid_dict[ftid]) > 1:
                # check if there is one of the multiple fusion_gene records where both partners correspond to the ftid
                for fusion_gene in ftid_dict[ftid]:
                    left_fusion, right_fusion = fusion_gene.split("_")
                    if left_fusion in ftid and right_fusion in ftid:
                        ftid_dict[ftid] = fusion_gene
                        break
                # if none of the fusion_gene records fit to the ftid, take the shortest and print a warning
                if type(ftid_dict[ftid]) == list:
                    print("Warning: Unresolvable annotation bias detected between ftid \"{0}\" and fusion_gene records \"{1}\"".format(ftid, ftid_dict[ftid]))
                    ftid_dict[ftid] = min(ftid_dict[ftid], key=len)
                    print("Warning: Chosing {} as fusion_gene record which might be wrong.".format(ftid_dict[ftid]))
            else:
                ftid_dict[ftid] = ftid_dict[ftid][0]

        # create sets of fusion gene names according to different fiters
        fusion_gene_set_list.append(set([ftid_dict[x] for x in df_rep[df_rep["type"] != "cis_near"]["FTID"]])) # exclude cis near events as they are likely read through
        fusion_gene_set_list.append(set([ftid_dict[x] for x in df_rep[df_rep["exon_boundary"] == "both"]["FTID"]])) # allow only fusions where the breakpoints are on exon boundaries in BOTH fusion partners
        fusion_gene_set_list.append(set([ftid_dict[x] for x in df_rep[df_rep["frame"] != "no_frame"]["FTID"]])) # exclude no frame fusions
        fusion_gene_set_list.append(set([ftid_dict[x] for x in df_rep[df_rep["neo_peptide_sequence"].str.isalpha()]["FTID"]])) # the neo peptide sequence must be an alphabetic sequence of at least length 1
        if is_merged:
            fusion_gene_set_list.append(set([ftid_dict[x] for x in df_rep[df_rep["ft_junc_cnt_org_1"] + df_rep["ft_span_cnt_org_1"] + df_rep["ft_junc_cnt_org_2"] + df_rep["ft_span_cnt_org_2"] > 0]["FTID"]])) # at least 1 junction or spanning counts on the fusion transcript in either replicate
        else:
            fusion_gene_set_list.append(set([ftid_dict[x] for x in df_rep[df_rep["ft_junc_cnt_org"] + df_rep["ft_span_cnt_org"] > 0]["FTID"]])) # at least 1 junction or spanning counts on the fusion transcript
        fusion_gene_set_list.append(set([ftid_dict[z] for z in df_rep[[all([not y.startswith(tuple(gene_blacklist)) for y in x]) for x in df_rep["Fusion_Gene"].str.split("_")]]["FTID"]])) # both fusion partners are not allowed to start with anything mentioned in the "blacklisted" list
        if is_merged:
            fusion_gene_set_list.append(set([ftid_dict[x] for x in df_rep[(df_rep["prediction_class_1"] == "positive") | (df_rep["prediction_class_2"] == "positive")]["FTID"]])) # the random forrest model must have classified this at least once as "posititve" in either replicate
        else:
            fusion_gene_set_list.append(set([ftid_dict[x] for x in df_rep[df_rep["prediction_class"] == "positive"]["FTID"]])) # the random forrest model must have classified this at least once as "posititve"
        fusion_gene_set_list.append(set([ftid_dict[x] for x in df_rep[(df_rep["exon_starts"] != "0") & (df_rep["exon_ends"] != "0")]["FTID"]])) # tmp hack to exclude valid ftids with missing exon numbers (bug in GetFusionSequence.R?)
        
        fusion_gene_filter_all = set(ftid_dict.values())
        fusion_gene_filter_wo_pred = set(ftid_dict.values())
        for i in range(len(fusion_gene_set_list)):
            fusion_gene_filter_all = fusion_gene_filter_all.intersection(fusion_gene_set_list[i])
            if i != 6:
                fusion_gene_filter_wo_pred = fusion_gene_filter_wo_pred.intersection(fusion_gene_set_list[i])
        
        fusion_gene_set_list.append(set(ftid_dict.values()))
        fusion_gene_set_list.append(fusion_gene_filter_all)
        fusion_gene_set_list.append(fusion_gene_filter_wo_pred)

        return (map(len, fusion_gene_set_list), fusion_gene_filter_all)
        # orginal implementation of the above which is much slower, less flexible, less pandas-style and more error prone
        # kept for testing only
#        for ftid in ftid_dict:
#            record2test = df_rep[df_rep["FTID"] == ftid]
#            # check whether one of the filtering criteria applies to this ftid
#            if(
#                    any(record2test["type"] == "cis_near") or
#                    any(record2test["exon_boundary"] != "both") or
#                    any(record2test["frame"] == "no_frame") or
#                    not all(record2test["neo_peptide_sequence"].str.isalpha()) or
#                    any(record2test["ft_junc_cnt_org"] + record2test["ft_span_cnt_org"] == 0) or
#                    any([x.startswith(tuple(gene_blacklist)) for x in record2test["Fusion_Gene"].values[0].split("_")]) or
#                    any(record2test["prediction_class"] == "negative") or
#                    any(record2test["exon_starts"] == "0") or
#                    any(record2test["exon_ends"] == "0") # see before
#                    ):
#                df_rep.drop(labels=record2test.index, inplace=True)
#        unfiltered_fg = set(ftid_dict.values())
#        filtered_fg = unfiltered_fg.intersection(set(df_rep["Fusion_Gene"]))

    @staticmethod
    def mimic_icam_implementation(inputfile1, inputfile2, outputfile):
        """ This is mimicing the fusion_peptide_filter.py implementation of MALO for testing purposes only! """
        fusion_filter_dict = {}
        id_filter_list = []
        selected_cols = ["FTID", "Fusion_Gene", "Breakpoint1", "Breakpoint2", "type", "exon_starts", "exon_ends", "exon_boundary",
                         "frame", "context_sequence", "neo_peptide_sequence", "neo_peptide_sequence_bp", "ft_junc_cnt_org", "ft_span_cnt_org",
                         "wt1_junc_cnt_org", "wt1_span_cnt_org", "wt2_junc_cnt_org", "wt2_span_cnt_org", "prediction_prob", "prediction_class",
                         "n_replicates"]
        cols_dict = dict(zip(selected_cols, range(0, 21)))
        
        for inputfile in [inputfile1, inputfile2]:
            # load required data from easyfuse output into pd dataframe and append column for replicate counting
            df_rep = pd.read_csv(inputfile, sep = ";")
            df_rep = df_rep[selected_cols[:-1]]
            #df_rep[str(selected_cols[-1])] = 0
            
            # iterate over rows in the df
            for i in df_rep.index:
                ftid = df_rep.loc[i, "FTID"]
                if ftid not in fusion_filter_dict:
                    fusion_filter_dict[ftid] = df_rep.loc[i].tolist() + [1]
                    if fusion_filter_dict[ftid][cols_dict["prediction_class"]].lower() == "positive":
                        fusion_filter_dict[ftid][cols_dict["prediction_class"]] = 1
                    else:
                        fusion_filter_dict[ftid][cols_dict["prediction_class"]] = 0
                    id_filter_list.append(ftid)
                else:
                    for target_read_cnt in ["ft_junc_cnt_org", "ft_span_cnt_org", "wt1_junc_cnt_org", "wt1_span_cnt_org", "wt2_junc_cnt_org", "wt2_span_cnt_org"]:
                        fusion_filter_dict[ftid][cols_dict[target_read_cnt]] += df_rep.loc[i, target_read_cnt]
                    if df_rep.loc[i, "prediction_class"].lower() == "positive":
                        fusion_filter_dict[ftid][cols_dict["prediction_class"]] += 1
                    fusion_filter_dict[ftid][cols_dict["prediction_prob"]] = max(fusion_filter_dict[ftid][cols_dict["prediction_prob"]], df_rep.loc[i, "prediction_prob"])
                    fusion_filter_dict[ftid][cols_dict["n_replicates"]] += 1
        # filter records from the dict
        id_passing_filter_list = []
        min_replicates = 2
        blacklisted = ["HLA", "IG"]
        for ftid in id_filter_list:
            if(
                    fusion_filter_dict[ftid][cols_dict["frame"]] != "no_frame" and # exclude no frame fusions
                    fusion_filter_dict[ftid][cols_dict["exon_boundary"]] == "both" and # allow only fusions where the breakpoints are on exon boundaries in BOTH fusion partners
                    fusion_filter_dict[ftid][cols_dict["prediction_class"]] >= 1 and # the random forrest model must have classified this at least once as "posititve"
                    fusion_filter_dict[ftid][cols_dict["ft_junc_cnt_org"]] + fusion_filter_dict[ftid][cols_dict["ft_span_cnt_org"]] > 0 and # at least 1 junction or spanning counts on the fusion transcript
                    all([not x.startswith(tuple(blacklisted)) for x in fusion_filter_dict[ftid][cols_dict["Fusion_Gene"]].split("_")]) and # both fusion partners are not allowed to start with everything mentioned in the "blacklisted" list
                    fusion_filter_dict[ftid][cols_dict["n_replicates"]] >= min_replicates and # number of input samples (i.e. files) with the same ftid (urla: this is not unique and should be ftid+context_sequence_100_id!)
                    fusion_filter_dict[ftid][cols_dict["neo_peptide_sequence"]].isalpha() and # the neo peptide sequence must be an alphabetic sequence of at least length 1
                    fusion_filter_dict[ftid][cols_dict["type"]] != "cis_near" and # exclude cis near events as they are likely read through
                    fusion_filter_dict[ftid][cols_dict["exon_starts"]] != "0" and # tmp hack to exclude valid ftids with missing exon numbers (bug in GetFusionSequence.R?)
                    fusion_filter_dict[ftid][cols_dict["exon_ends"]] != "0" # see before
                    ):
                id_passing_filter_list.append(ftid)
        # write filtered records to file
        with open(outputfile, "w") as outfile:
            outfile.write("{}\n".format("\t".join(selected_cols)))
            for ftid in id_passing_filter_list:
                outfile.write("{}\n".format("\t".join([str(x) for x in fusion_filter_dict[ftid]])))

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
