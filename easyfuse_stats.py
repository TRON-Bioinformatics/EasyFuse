#!/usr/bin/env python
"""
TODO: doc
"""
from __future__ import print_function
import re
import os
import os.path
import sys
import json
from datetime import datetime
import time
#import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import seaborn as sns
from easyfuse_samples import Samples
from easyfuse_scheduling import Queue
from easyfuse_statscollect import SCollector

#class Stats():
#def __init__(self):
#input_path=/mnt/bfx/urla/PRJNA252360_syntheticSpikeInFusions/fusionProcessing

class Stats(object):
    """Collect stats of the run and write them to file"""
    # pylint: disable=too-many-instance-attributes
    #         min_rl_perc might be removed, other are required
    def __init__(self, input_path, is_cpu_test):
        self.input_path = input_path
        self.is_cpu_test = is_cpu_test
        self.sample = Samples(os.path.join(input_path, "samples.csv"))
        self.sample_stats = SCollector(os.path.join(input_path, "sample_stats.tdt"))
        self.sample_read_cnts = {}
        self.sample_read_cnts_org = {} # original input read counts
        self.sample_read_cnts_flt = {} # filtered read counts
#        self.max_read_cnt = 0
        self.figure_list = []

    def plot_filter_stats(self):
        """Parses read filtering log file in a dictionary for plotting"""
        filter_dict_raw = {} # dict of raw/absolute values
        filter_dict_all_pc = {} # dict of relative values (percentages); only MM, fltr and other ("normal") plus QC1/2
        filter_dict_fltr_pc = {} # dict of relative values (percentages); distribution of all filtered groups (chim, un, npp, clip) plus QC1/2
        for sample_id in self.sample.sample_map:
            if not "Readfilter" in self.sample.get_tool_list_from_state(sample_id):
                filter_dict_raw[sample_id] = [-1, -1, -1, -1]
                filter_dict_all_pc[sample_id] = [-1, -1, -1, -1]
                filter_dict_fltr_pc[sample_id] = [-1, -1, -1, -1]
                continue
            filtered_pairs = -1
            allinput_pairs = -1
            chimeric_pairs = -1
            multimap_pairs = -1
            unmapped_pairs = -1
            discorda_pairs = -1
            clipping_pairs = -1
            nochimer_pairs = -1
            qc1_test = -1
            qc2_test = -1
            zqc_test = -1
            file_open = os.path.join(self.input_path, "Sample_{}".format(sample_id), "scratch", "filtered_reads", "{}_Aligned.out.filtered.bam.fusionReadFilterLog".format(sample_id))
            with open(file_open, "r") as filter_log:
                for line in filter_log:
                    if "INFO Processed" in line:
                        pattern_match = re.search(r'alignments, ([0-9]*) of ([0-9]*) pairs', line)
                        filtered_pairs = float(pattern_match.group(1))
                        allinput_pairs = float(pattern_match.group(2))
                    elif "INFO Star_chimeric" in line:
                        chimeric_pairs = float(re.search(r'(\t|: )([0-9]*) pairs ', line).group(2))
                    elif "INFO Multimapped" in line:
                        multimap_pairs = float(re.search(r'(\t|: )([0-9]*) pairs ', line).group(2))
                    elif "INFO Unmapped" in line:
                        unmapped_pairs = float(re.search(r'(\t|: )([0-9]*) pairs ', line).group(2))
                    elif "INFO No proper pair" in line:
                        discorda_pairs = float(re.search(r'(\t|: )([0-9]*) pairs ', line).group(2))
                    elif "INFO 10bp_s_clip" in line or "INFO 10bp_non_M_clip" in line: # last part is for compatibility with older Readfilter outputs
                        clipping_pairs = float(re.search(r'(\t|: )([0-9]*) pairs ', line).group(2))
                    elif "INFO Unlikely chimeric" in line:
                        nochimer_pairs = float(re.search(r'(\t|: )([0-9]*) pairs ', line).group(2))
                    elif "INFO Filter QC1" in line:
                        pattern_match = re.search(r'(\t|: )(True|False)', line).group(2)
                        if pattern_match == "True":
                            qc1_test = 1
                    elif "INFO Filter QC2" in line:
                        pattern_match = re.search(r'(\t|: )(True|False)', line).group(2)
                        if pattern_match == "True":
                            qc2_test = 1
            # combined qc1/2 expression. if only qc1 is true, zqc becomes 0.66, if only qc2, zqc becomes 0.33; if both 1, if none 0
            if qc1_test == 1:
                if qc2_test == 1:
                    zqc_test = 1
                else:
                    zqc_test = 0.66
            else:
                if qc2_test == 1:
                    zqc_test = 0.33
                else:
                    zqc_test = 0

            self.sample_read_cnts_org[sample_id] = allinput_pairs
            self.sample_read_cnts_flt[sample_id] = filtered_pairs

            filter_dict_raw[sample_id] = [filtered_pairs, allinput_pairs/100, chimeric_pairs, multimap_pairs, unmapped_pairs, discorda_pairs, clipping_pairs]
            filter_dict_all_pc[sample_id] = [filtered_pairs / allinput_pairs, multimap_pairs / allinput_pairs, nochimer_pairs / allinput_pairs, zqc_test]
            filter_dict_all_pc[sample_id] = [round(100 * x, 2) for x in filter_dict_all_pc[sample_id]]
            filter_dict_fltr_pc[sample_id] = [chimeric_pairs / filtered_pairs, unmapped_pairs / filtered_pairs, discorda_pairs / filtered_pairs, clipping_pairs / filtered_pairs, zqc_test]
            filter_dict_fltr_pc[sample_id] = [round(100 * x, 2) for x in filter_dict_fltr_pc[sample_id]]
        self.plot_heatmap(filter_dict_raw, ["Fltr", "In(x100)", "Chim", "MM", "UN", "NPP", "Clip"], ["Sample", "PairType", "Data"], "Filter counts of read filtering")
        self.plot_heatmap(filter_dict_all_pc, ["Fltr", "MM", "Norm", "ZQC"], ["Sample", "PairType", "Data"], "Read filtering summary")
        self.plot_heatmap(filter_dict_fltr_pc, ["Chim", "UN", "NPP", "Clip", "ZQC"], ["Sample", "PairType", "Data"], "Filtered pairs summary")

    def plot_star_stats(self):
        """Parses a star output log file in a dictionary for plotting"""
        star_dict = {}
        # for each sample, get the data from the star_Log.final.out file
        # Store them in a dict and send it to plot_stacked_bars for plotting
        for sample_id in self.sample.sample_map:
            # set all values to -1 if star hasn't run (successfully)
            if not "Star" in self.sample.get_tool_list_from_state(sample_id):
                star_dict[sample_id] = [-1, -1, -1, -1]
                continue
            uniquely_mapped = 0
            multi_mapped = 0
            unmapped = 0
            chimeric = 0
            file_open = os.path.join(self.input_path, "Sample_{}".format(sample_id), "scratch", "expression", "star", "{}_Log.final.out".format(sample_id))
            with open(file_open, "r") as star_log:
                for line in star_log:
                    left_match = line.split("|")[0].strip()
                    if left_match == "Number of input reads" and not sample_id in self.sample_read_cnts:
                        self.sample_read_cnts[sample_id] = line.split("|")[1].strip()
                    if left_match == "Uniquely mapped reads %":
                        uniquely_mapped = float(line.split("|")[1].strip()[:-1])
                    elif left_match == "% of reads mapped to multiple loci" or left_match == "% of reads mapped to too many loci":
                        multi_mapped += float(line.split("|")[1].strip()[:-1])
                    elif left_match == "% of reads unmapped: too many mismatches" or left_match == "% of reads unmapped: too short" or left_match == "% of reads unmapped: other":
                        unmapped += float(line.split("|")[1].strip()[:-1])
                    elif left_match == "% of chimeric reads":
                        chimeric = float(line.split("|")[1].strip()[:-1])
            # add values to the data dict
            star_dict[sample_id] = [uniquely_mapped, multi_mapped, unmapped, chimeric]
            # allow max two decimal points 
            star_dict[sample_id] = [round(x, 2) for x in star_dict[sample_id]]
#        self.plot_stacked_bars(star_dict, ["UMR", "MMR", "UNR", "CMR"], ("Uniquely-, multi-, un- and chimeric mapped reads by STAR", "", "Read mappings (%)"), True, True)
        self.plot_heatmap(star_dict, ["UMR", "MMR", "UNR", "CMR"], ["Sample", "Mapping", "Data"], "Uniquely-, multi-, un- and chimeric mapped reads by STAR")

    def plot_kallisto_stats(self):
        """Parses a kallisto output file in a dictionary for plotting"""
        kallisto_dict = {}
        # for each sample, get the data from the run_info.json file
        # Store them in a dict and send it to plot_stacked_bars for plotting
        for sample_id in self.sample.sample_map:
            # set all values to -1 if kallisto hasn't run (successfully)
            if not "Kallisto" in self.sample.get_tool_list_from_state(sample_id):
                kallisto_dict[sample_id] = [-1, -1, -1]
                continue
            file_open = os.path.join(self.input_path, "Sample_{}".format(sample_id), "scratch", "expression", "kallisto", "run_info.json")
            with open(file_open, "r") as info_file:
                json_data = json.load(info_file)
            uniquely_mapped = json_data["p_unique"]
            multi_mapped = json_data["p_pseudoaligned"] - json_data["p_unique"]
            unmapped = 100 - json_data["p_pseudoaligned"]
            kallisto_dict[sample_id] = [uniquely_mapped, multi_mapped, unmapped]
            # allow max two decimal points 
            kallisto_dict[sample_id] = [round(x, 2) for x in kallisto_dict[sample_id]]
            
            # get read counts and identify the max read cnts for normalizations afterwards
            processed_reads = json_data["n_processed"]
            if not sample_id in self.sample_read_cnts:
                self.sample_read_cnts[sample_id] = processed_reads
#            self.max_read_cnt = max(self.max_read_cnt, processed_reads)
#        self.plot_stacked_bars(kallisto_dict, ["UMR", "MMR", "UNR"], ("Uniquely-, multi- and unmapped reads by Kallisto", "", "Read mappings (%)"), True, True)
        self.plot_heatmap(kallisto_dict, ["UMR", "MMR", "UNR"], ["Sample", "Mapping", "Data"], "Uniquely-, multi- and unmapped reads by Kallisto")

    def plot_fusion_number_stats(self):
        """Parses a kallisto output file in a dictionary for plotting"""
        fusion_number_dict = {}
        fusion_tool_list = ["Fusioncatcher", "Mapsplice", "Pizzly", "Starchip", "Starfusion", "Infusion"]
        # for each sample, get the data from the run_info.json file
        # Store them in a dict and send it to plot_stacked_bars for plotting
        for sample_id in self.sample.sample_map:
            fusion_root_path = os.path.join(self.input_path, "Sample_{}".format(sample_id), "scratch", "fusion")
            fusioncatcher_count = self.get_line_number_from_file(os.path.join(fusion_root_path, "fusioncatcher", "final-list_candidate-fusion-genes.txt"), True)
            mapsplice_count = self.get_line_number_from_file(os.path.join(fusion_root_path, "mapsplice", "fusions_well_annotated.txt"), False)
            pizzly_count = self.get_line_number_from_file(os.path.join(fusion_root_path, "pizzly", "kallizzy.json.txt"), True)
            starchip_count = self.get_line_number_from_file(os.path.join(fusion_root_path, "starchip", "starchip.summary"), True)
            starfusion_count = self.get_line_number_from_file(os.path.join(fusion_root_path, "starfusion", "star-fusion.fusion_predictions.abridged.tsv"), True)
            infusion_count = self.get_line_number_from_file(os.path.join(fusion_root_path, "infusion", "fusions.detailed.txt"), True)
            fusion_number_dict[sample_id] = [fusioncatcher_count, mapsplice_count, pizzly_count, starchip_count, starfusion_count, infusion_count]
            # write counts to stats file
            for i, fusion_count in enumerate(fusion_number_dict[sample_id]):
                self.sample_stats.add_stat(sample_id, fusion_tool_list[i], "fusion_count", fusion_count)

        # send data to plotting methods
#        self.plot_stacked_bars(fusion_number_dict, fusion_tool_list, ("Identified gene fusions per tool","","fusion events (#)"), False, True)
        self.plot_heatmap(fusion_number_dict, fusion_tool_list, ["Sample", "FusionTool", "Data"], "Identified gene fusions per tool")
        self.plot_scatter_correlation(self.sample_read_cnts, fusion_number_dict, fusion_tool_list, ['sampleID', 'fusiontool', 'reads', 'fusCount'], ("", "", ""))
    
    def plot_fusion_genes_data(self):
        for sample_id in self.sample.sample_map:
            fetched_fusions_path = os.path.join(self.input_path, "Sample_{}".format(sample_id), "fetchdata", "fd_1_tool", "fetched_fusions")
            if os.path.isdir(fetched_fusions_path):
                fusion_list_dict = {}
                fusion_files = [f for f in os.listdir(fetched_fusions_path) if os.path.isfile(os.path.join(fetched_fusions_path, f)) and "_toFusionInspector.txt" in f]
                for data_file in fusion_files:
                    gene_fusion_list = []
                    gene_fusion_count = 0
                    with open(os.path.join(fetched_fusions_path, data_file), "r") as fufi:
                        for gene_fusion_count, line in enumerate(fufi, 1):
                            gene_fusion_list.append(line.rstrip("\n"))
                    
                    fusion_list_dict[data_file.split("_")[0]] = gene_fusion_list
                    self.sample_stats.add_stat(sample_id, data_file.split("_")[0], "fusion_gene_count", "all:{}, uniq:{}".format(gene_fusion_count, len(set(gene_fusion_list))))
                    self.sample_stats.add_stat(sample_id, data_file.split("_")[0], "fusion_gene_list", set(gene_fusion_list))
                
                # compare results (all vs all) overlap from the fusion tools                
                for fusion_tool in fusion_list_dict:
                    comparison_list = []
                    for tool1_fusion in fusion_list_dict[fusion_tool]:
                        comparison_list.append(0)
                    for fusion_tool_vs in fusion_list_dict:
                        if fusion_tool == fusion_tool_vs:
                            continue
                        # compare predictions
                        count_hits = 0
                        for fus_pos, tool1_fusion in enumerate(fusion_list_dict[fusion_tool]):
                            tool1_geneA, tool1_geneB = tool1_fusion.split("--")
                            for tool2_fusion in fusion_list_dict[fusion_tool_vs]:
                                tool2_geneA, tool2_geneB = tool2_fusion.split("--")
                                if (tool1_geneA == tool2_geneA and tool1_geneB == tool2_geneB) or (tool1_geneB == tool2_geneA and tool1_geneA == tool2_geneB):
                                #if tool1_fusion in fusion_list_dict[fusion_tool_vs]:
                                    count_hits += 1
                                    comparison_list[fus_pos] += 1
                                
                        self.sample_stats.add_stat(
                                sample_id, 
                                fusion_tool, 
                                "fusion_overlap_with_{}".format(fusion_tool_vs), 
                                "{0:.2f}".format(float(count_hits) / float(len(fusion_list_dict[fusion_tool])))
                                )
                    for i in range(max(comparison_list) + 1):
                        self.sample_stats.add_stat(
                                sample_id, 
                                fusion_tool, 
                                "fusion_overlap_with_{}_tools".format(i), 
                                "{0:.2f}".format(float(comparison_list.count(i)) / float(len(fusion_list_dict[fusion_tool])))
                                )
                        
            #self.get_fusion_genes_data_file_list(sample_id, fetched_fusions_path)

    def get_fusion_genes_data_file_list(self, sample_id, fetched_fusions_path):
        if os.path.isdir(fetched_fusions_path):
            fusion_files = [f for f in os.listdir(fetched_fusions_path) if os.path.isfile(os.path.join(fetched_fusions_path, f)) and "_toFusionInspector.txt" in f]
            for data_file in fusion_files:
                gene_fusion_list = []
                gene_fusion_count = 0
                with open(os.path.join(fetched_fusions_path, data_file), "r") as fufi:
                    for gene_fusion_count, line in enumerate(fufi):
                        gene_fusion_list.append(line)
                        
                self.sample_stats.add_stat(sample_id, data_file.split("_")[0], "fusion_gene_count", gene_fusion_count)

    def plot_run_time_stats(self, tool_list_to_check, sact_out_dict):
        """Gets run times of slurm scripts via sacct and plots them"""
        run_time_dict = {}
        tool_list_to_plot = []
        for i, sample_id in enumerate(sact_out_dict):
            sample_run_time_list = []
            # walk through the ordered list of tools to ensure correct assignment of values
            # because tools can be differently sorted in the sample.csv file and in the sacct output
            for tool in tool_list_to_check:
                found_tool = 0
                for sacct_line in sact_out_dict[sample_id]:
                    if sacct_line.strip().startswith("{}-".format(tool)): # not adding the hyphen in the comparison would match for example star also to starfusion
                    #if sacct_line and not sacct_line.strip().startswith("JobName") and not sacct_line.strip().startswith("-"):
                        found_tool += 1
                        time_splitter = sacct_line.split()[1].split(":")
                        time_splitter = [int(j) for j in time_splitter]
                        tool_time = time_splitter[0] * 3600 + time_splitter[1] * 60 + time_splitter[2]
                        self.sample_stats.add_stat(sample_id, tool, "runtime", tool_time)
                        sample_run_time_list.append(tool_time)
                        if i == 0:
                            if self.is_cpu_test:
                                tool_list_to_plot.append("{}_{}".format(tool, sacct_line.split()[0].split("-")[2]))
                            else:
                                tool_list_to_plot.append(tool)
                # check the stats file if it wasn't found in the sacct data (maybe already deleted?)
                if found_tool == 0:
                    #self.sample_stats.
                    sample_run_time_list.append(0)
                    if i == 0:
                        tool_list_to_plot.append(tool)
            run_time_dict[sample_id] = sample_run_time_list
        # plot absolute values
#        self.plot_stacked_bars(run_time_dict, tool_list_to_check, ("Runtime of individual programs","","runtime (s)"), False, False)
        self.plot_heatmap(run_time_dict, tool_list_to_plot, ["Sample", "Program", "Data"], "Runtime of inidividual programs")
        self.plot_scatter_correlation(self.sample_read_cnts, run_time_dict, tool_list_to_plot, ['sampleID', 'fusiontool', 'reads', 'run_time'], ("", "", ""))
        print("prev read cnts: {}".format(self.sample_read_cnts_org))
        print("new org read cnts: {}".format(self.sample_read_cnts_org))
        print("new flt read cnts: {}".format(self.sample_read_cnts_flt))

    def plot_heatmap(self, data_dict, group_names, column_names, chart_title):
        """Create a pandas dataframe from a dictionary and plot the data table as a heatmap.
        Summarize the data furthermore in a boxplot with swarmplot overlay"""
#        print(data_dict)
#        print(group_names)
        list_of_data_lists = []
        min_data_value = sys.maxsize
        max_data_value = 0
        for data_dict_key in data_dict:
            for i, list_record in enumerate(data_dict[data_dict_key]):
                min_data_value = min(list_record, min_data_value)
                max_data_value = max(list_record, max_data_value)
                list_of_data_lists.append([data_dict_key, group_names[i], list_record])
        if len(data_dict) > 20:
            sub_figure_list = range(0, len(data_dict)*len(group_names), 20*len(group_names))
            sub_figure_list.append(len(data_dict)*len(group_names))
            for i in range(0, len(sub_figure_list) - 1):
                sub_list_of_data_lists = list_of_data_lists[sub_figure_list[i]:sub_figure_list[i+1]]
                sub_list_of_data_lists.append(["MinMax", group_names[0], min_data_value])
                sub_list_of_data_lists.append(["MinMax", group_names[1], max_data_value])
                for j in range(2, len(group_names)):
                    sub_list_of_data_lists.append(["MinMax", group_names[j], 0])
                self.plot_heatmap_helper(sub_list_of_data_lists, column_names, chart_title)
        else:
            self.plot_heatmap_helper(list_of_data_lists, column_names, chart_title)
        # plot data additionally as box plot (plus swarmplot overlay)
        pd_df = pd.DataFrame(list_of_data_lists, columns=column_names)
        plt.figure(len(self.figure_list) + 1)
        ax = plt.axes()
        plt.title(chart_title)
        sns.boxplot(data=pd_df, ax=ax, x=column_names[1], y=column_names[2], hue=column_names[1], dodge=False)
        sns.swarmplot(data=pd_df, ax=ax, x=column_names[1], y=column_names[2], color=(0.8, 0.8, 0.8), alpha=0.5)
        plt.tight_layout()
        self.figure_list.append(ax.figure)

    def plot_heatmap_helper(self, list_of_data_lists, column_names, chart_title):
        """actual dataframe generation and heatmap plotting"""
        # create a pandas data frame with one column for each sample names, program(oe) names and data values
        pd_df = pd.DataFrame(list_of_data_lists, columns=column_names)
        pd_df_pivot = pd_df.pivot(index=column_names[0], columns=column_names[1], values=column_names[2])
        plt.figure(len(self.figure_list) + 1)
        ax = plt.axes()
        plt.title(chart_title)
        sns.heatmap(pd_df_pivot, ax=ax, annot=True, fmt="g", annot_kws={"size": 4}, cmap="tab20b", yticklabels=True)
        plt.tight_layout()
        self.figure_list.append(ax.figure)

    def plot_scatter_correlation(self, x_data_dict, y_data_dict, data_aggregate_list, column_names, chart_titles):
        """Combine two corresponding data dicts into a pd data frame and create+save a scatter plot from it"""
        # column_names are expected to be "sample_id" / "data_aggregation_type" (e.g. fusiontool) / "data for x axis" / "data for y axis"
        final_df = pd.DataFrame(columns=column_names)
        counter = 0
        for i, fus_tool in enumerate(data_aggregate_list, 0):
            for sample_id in self.sample.sample_map:
                if sample_id in y_data_dict and sample_id in x_data_dict:
                    # fast hack to correct runtime comparison for "Readfilter"
                    if fus_tool == "Readfilter" and "run_time" in data_aggregate_list:
                        final_df.loc[counter] = [sample_id, fus_tool, self.sample_read_cnts_flt[sample_id], y_data_dict[sample_id][i]]
                    else:
                        final_df.loc[counter] = [sample_id, fus_tool, x_data_dict[sample_id], y_data_dict[sample_id][i]]
                    counter += 1
        final_df[[column_names[2], column_names[3]]] = final_df[[column_names[2], column_names[3]]].apply(pd.to_numeric)
        sns_lmplot = sns.lmplot(x=column_names[3], y=column_names[2], col=column_names[1], hue=column_names[0], data=final_df, col_wrap=2, sharex=False, sharey=False, fit_reg=False, scatter_kws={"s": 100, "alpha": 1}, palette="muted")
#        sns_lmplot.fig.savefig(save_path)
        self.figure_list.append(sns_lmplot.fig)

#    def list_predicted_gene_fusions(self, fusion_tool_cutoff):
#        #Sample_SRR1659951/fetchdata/fd_2_tool/fetched_fusions> head Detected_Fusions.csv
#        fusion_dict = {}
#        for sample_id in self.sample.sample_map:
#            if not "Fetchdata" in self.sample.get_tool_list_from_state(sample_id):
#                print("It looks like fusions weren't collected - did fetchdata finish properly?! Skipping {}".format(sample_id))
#                continue
#            detect_fus_file = os.path.join(self.input_path, "Sample_{}".format(sample_id), fetchdata, "fd_{}_tool".format(fusion_tool_cutoff), "fetched_fusions", "Detected_Fusions.csv")
#            with open(detect_fus_file, "r") as csv_file:                            
            
    @staticmethod
    def get_line_number_from_file(in_file, has_header):
        """Return the number of lines of a file"""
        try:
            i = 0
            with open(in_file) as test_file:
                for i, _ in enumerate(test_file):
                    pass
            if has_header:
                return i
            return i + 1
        except IOError:
            print("File {} is not present or not accessible".format(in_file))
            return -1

    def get_sample_processing_info(self, save_path):
        """Write a list of finished slurm jobs sorted by their submission date"""
        job_id_dict = {}
        for sample_id in self.sample.sample_map:
            if not "," in self.sample.get_state_string(sample_id):
                continue
            for field in self.sample.get_state_string(sample_id).split(", "):
                field_splitter = field.split("-")
                job_id_dict[field_splitter[2]] = (field_splitter[1], field_splitter[0])
        #OrderedDict(sorted(job_id_dict.items(), key=lambda t: t[0]))
        with open(save_path, "w") as out_file:
            last_day = ""
            for sort_samp_id in sorted(job_id_dict.keys()):
                current_day = datetime.fromtimestamp(float(sort_samp_id)).strftime("%Y-%m-%d")
                if not current_day == last_day:
                    out_file.write("{}\n".format(current_day))
                    last_day = current_day        
                out_file.write("{}\t{}\t{}\n".format(datetime.fromtimestamp(float(sort_samp_id)).strftime("%H:%M:%S"), job_id_dict[sort_samp_id][0], job_id_dict[sort_samp_id][1]))

    def run(self):
        """Execute individual methods"""
        if not self.is_cpu_test:
            print("Plotting filtering stats...")
            self.plot_filter_stats()
            print("Plotting star mapping stats...")
            self.plot_star_stats()
            print("Plotting kallisto mapping stats...")
#            self.plot_kallisto_stats()
            print("Plotting fusion number stats...")
            self.plot_fusion_number_stats()
            print("Plotting fusion number stats2...")
            self.plot_fusion_genes_data()
        print("Trying to get slurm stats for all samples...")
        sact_out_dict = {}
        for sample_id in self.sample.sample_map:
            # get the complete state string of a sample for querying sacct and 
            # (derived from the same str), the run tools as a list
            max_state_len = 0
            sample_state = self.sample.get_state_string(sample_id)
            for split_state in sample_state.split(", "):
                max_state_len = max(max_state_len, len(split_state))
            # get job information with sacct
            cmd_timeCollect = "sacct --name {0} --allocations --format=JobName%{1},Elapsed,Start,ExitCode,State,AveCPU,AveVMSize,MaxVMSize,MaxRSS,reqmem,NCPUS,Priority,ReqCPUS,Reserved".format(sample_state.replace(", ", ","), max_state_len)
            sacct_stdout_string = Queue.submit_nonqueue_and_get_stdout(cmd_timeCollect.split(" "))[0]
            sact_out_dict[sample_id] = sacct_stdout_string.split("\n")
        print("Plotting run time stats...")
        self.plot_run_time_stats(["Readfilter", "Kallisto", "Star", "Fusioncatcher", "Mapsplice", "Pizzly", "Starfusion", "Starchip", "Infusion"], sact_out_dict)
        now_time = datetime.fromtimestamp(time.time()).strftime("%Y%m%d_%H%M%S")
        self.get_sample_processing_info(os.path.join(self.input_path, "RunStats_{}_slurm_submissions_of_finished_jobs.txt".format(now_time)))

        with PdfPages(os.path.join(self.input_path, "RunStats_{}.pdf".format(now_time))) as pdf:
            for plot_fig in self.figure_list:
                pdf.savefig(plot_fig)

def main():
    """Command line argument parsing and app start"""
    parser = argparse.ArgumentParser(description='Post processing of an easyfuse run - currently, collecting runtimes only :)')

    parser.add_argument('-i', '--input', dest='input', help='Specify the easyfuse root dir of the run you want to process.', required=True)
    parser.add_argument('--cpu_test', dest='cpu_test', help=argparse.SUPPRESS, default=False, action='store_true')
    args = parser.parse_args()

    stats = Stats(args.input, args.cpu_test)
    stats.run()

if __name__ == '__main__':
    main()
