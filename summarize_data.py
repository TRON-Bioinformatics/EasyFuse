#!/usr/bin/env python

"""
@author: BNT (URLA)
@version: 20190118
"""

from __future__ import print_function
import os
import os.path
from datetime import datetime
import time
#import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import seaborn as sns
from misc.config import Config
from join_data import DataJoining
from misc.samples import Samples
import misc.io_methods as IOMethods

class IcamSummary(object):
    """Collect stats of the run and write them to file"""
    def __init__(self, input_path, config_path):
        self.input_path = input_path
        self.samples = Samples(os.path.join(input_path, "samples.db"))
        self.cfg = Config(config_path)
        self.figure_list = []

    def run(self, icam_run, model_predictions):
        """Execute individual methods"""
        fusion_tools = self.cfg.get("general", "fusiontools").split(",")
        fusion_data_summary_path = os.path.join(self.input_path, "FusionSummary")
        IOMethods.create_folder(fusion_data_summary_path)

        now_time = datetime.fromtimestamp(time.time()).strftime("%Y%m%d_%H%M%S")
        figure_out_path = os.path.join(fusion_data_summary_path, "FusionSelection_figures_{}.pdf".format(now_time))
        data_out_path_fltr = os.path.join(fusion_data_summary_path, "FusionSelection_data_filtering_{}.txt".format(now_time))
        data_out_path_freq_all = os.path.join(fusion_data_summary_path, "FusionSelection_data_fusFreqsSingles_{}.txt".format(now_time))
        data_out_path_freq_12 = os.path.join(fusion_data_summary_path, "FusionSelection_data_fusFreqsMerged_{}.txt".format(now_time))

        fusion_frequency_all = {}
        fusion_frequency_12 = {}
        filtering_data_1 = {}
        filtering_data_2 = {}
        filtering_data_12 = {}
        # 0:type; 1:boundary; 2:frame; 3:pepseq; 4:counts; 5:blacklist; 6:prediction;
        # 7:exonError; 8:unfiltered; 9:allFilter; 10:allButPredFilter
        colnames_filtering = [
                "cis_near_no", "exon_bound_both", "no_frame_no", "pep_seq_yes", "counts_yes", "blacklist_no", "pred_pos",
                "exon_err_no", "unfiltered", "all_fltr", "all_fltr_but_pred"
                ]
        # get sorted list of sample ids
        sid_list = sorted(self.samples.get_sample_id_list())
#        sid_list = sid_list[:10]
        current_base = ""
        sample1 = ""
        sample2 = ""

        count_pairs = 0
        count_valid_pairs = 0
        count_valid_sample = 0
        # check how many samples are available for processing and create {sample_id: sample_date} dict
        sample_date_dict = {}
        sample_toolCnt_dict = {}
        for i, sample in enumerate(sid_list, 1):
            try:
                sample_date_dict[sample] = self.samples.get_fastq_files(sample)[0].split("/")[8].split("_")[0]
            except IndexError:
                print("Warning: Could not retrieve sample date for {}".format(sample))
                sample_date_dict[sample] = "NA"
            # count the numbers of fusion tools that have been run on this sample
            sample_toolCnt_dict[sample] = len([tool for tool in fusion_tools if tool in self.samples.get_tool_list_from_state(sample)])
            if "Fetchdata" in self.samples.get_tool_list_from_state(sample):
                count_valid_sample += 1
            if icam_run:
                if "Rep1" in sample:
                    current_base = sample.split("_")[0]
                    sample1 = sample
                if "Rep2" in sample and current_base == sample.split("_")[0]:
                    sample2 = sample
                    count_pairs += 1
                    if sample_date_dict[sample1] == sample_date_dict[sample2]:
                        sample_date_dict[current_base] = sample_date_dict[sample1]
                    else:
                        sample_date_dict[current_base] = sample_date_dict[sample1] + "_" + sample_date_dict[sample2]
                    # count the numbers of fusion tools that have been run on either sample1 or sample2
                    current_base_fus_set = set(map(str, self.samples.get_tool_list_from_state(sample1)) + map(str, self.samples.get_tool_list_from_state(sample2)))
                    sample_toolCnt_dict[current_base] = len([tool for tool in fusion_tools if tool in current_base_fus_set])
                    if "Fetchdata" in self.samples.get_tool_list_from_state(sample1) and "Fetchdata" in self.samples.get_tool_list_from_state(sample2):
                        count_valid_pairs += 1
        if icam_run:
            print("Found {0} (partially) processed samples (from {1} patients) in {2}. Data will be collected from {3} patient samples for which fetchdata has been run on both replicates.".format(i, count_pairs, self.input_path, count_valid_pairs))
        else:
            print("Found {0} (partially) processed samples in {1}. Data will be collected from {2} samples for which fetchdata has been run.".format(i, self.input_path, count_valid_sample))

        average_time = 0
        count_processed = 0
        with open(data_out_path_fltr, "w") as fltr_out:
            fltr_out.write("SampleID\tSampleDate\tFusionToolCnt\t" + "\t".join(colnames_filtering) + "\n")
            for sample in sid_list:
                if icam_run:
                    if "Rep1" in sample:
                    #if "S3" in sample or "S5" in sample:
                        current_base = sample.split("_")[0]
                        sample1 = sample
                    if "Rep2" in sample and current_base == sample.split("_")[0]:
                    #if ("S4" in sample or "S6" in sample) and current_base == sample.split("_")[0]:
                        sample2 = sample
                        #call fusionranking
                        if "Fetchdata" in self.samples.get_tool_list_from_state(sample1) and "Fetchdata" in self.samples.get_tool_list_from_state(sample2):
                            count_processed += 1
                            print("Processing patient ID {0} (dataset {1}/{2})".format(current_base, count_processed, count_valid_pairs))
                            start_time = time.time()
                            fusion_data_summary = DataJoining(self.input_path, sample1, sample2, os.path.join(fusion_data_summary_path, current_base), False).run(self.cfg, icam_run, model_predictions)

                            fusion_frequency_all = self.add_to_fus_dict(fusion_data_summary[0][1], fusion_frequency_all)
                            fusion_frequency_all = self.add_to_fus_dict(fusion_data_summary[1][1], fusion_frequency_all)
                            fusion_frequency_12 = self.add_to_fus_dict(fusion_data_summary[2][1], fusion_frequency_12)
                            
                            filtering_data_1[sample1] = fusion_data_summary[0][0]
                            filtering_data_2[sample2] = fusion_data_summary[1][0]
                            filtering_data_12[current_base] = fusion_data_summary[2][0]
                            
                            # append filtering plots per sample
                            self.plot_boxswarm({sample1:fusion_data_summary[0][0]}, colnames_filtering, ["Sample", "Filter", "Data"], "Filter counts of distinct fusion genes for sample {} - Rep1's".format(sample1))
                            self.plot_boxswarm({sample2:fusion_data_summary[1][0]}, colnames_filtering, ["Sample", "Filter", "Data"], "Filter counts of distinct fusion genes for sample {} - Rep2's".format(sample2))
                            self.plot_boxswarm({current_base:fusion_data_summary[2][0]}, colnames_filtering, ["Sample", "Filter", "Data"], "Filter counts of distinct fusion genes for sample {} - Merged".format(current_base))

                            fltr_out.write(sample1 + "\t" + str(sample_date_dict[sample1]) + "\t" + str(sample_toolCnt_dict[sample1]) + "\t" + "\t".join(map(str, fusion_data_summary[0][0])) + "\n")
                            fltr_out.write(sample2 + "\t" + str(sample_date_dict[sample2]) + "\t" + str(sample_toolCnt_dict[sample2]) + "\t" + "\t".join(map(str, fusion_data_summary[1][0])) + "\n")
                            fltr_out.write(current_base + "\t" + str(sample_date_dict[current_base]) + "\t" + str(sample_toolCnt_dict[current_base]) + "\t" + "\t".join(map(str, fusion_data_summary[2][0])) + "\n")

                            time_taken = time.time() - start_time
                            average_time = (average_time * (count_processed-1) + time_taken) / count_processed
                            estimated_end = average_time * (count_valid_pairs - count_processed)
                            print("done. Processing time: {0:.2f}s; Average processing time: {1:.2f}s; Estimated end of processing in {2:.2f}s)".format(time_taken, average_time, estimated_end))
                else:
                    if "Fetchdata" in self.samples.get_tool_list_from_state(sample):
                        count_processed += 1
                        print("Processing sample {0} (dataset {1}/{2})".format(sample, count_processed, len(sid_list)))
                        start_time = time.time()
                        fusion_data_summary = DataJoining(self.input_path, sample, "", os.path.join(fusion_data_summary_path, sample), False).run(self.cfg, icam_run, model_predictions)

                        fusion_frequency_all = self.add_to_fus_dict(fusion_data_summary[1], fusion_frequency_all)
                        filtering_data_1[sample] = fusion_data_summary[0]
                        self.plot_boxswarm({sample:fusion_data_summary[0]}, colnames_filtering, ["Sample", "Filter", "Data"], "Filter counts of distinct fusion genes for sample {}".format(sample))
                        
                        fltr_out.write(sample + "\t" + str(sample_date_dict[sample]) + "\t" + str(sample_toolCnt_dict[sample]) + "\t" + "\t".join(map(str, fusion_data_summary[0])) + "\n")

                        time_taken = time.time() - start_time
                        average_time = (average_time * (count_processed-1) + time_taken) / count_processed
                        estimated_end = average_time * (count_valid_sample - count_processed)
                        print("done. Processing time: {0:.2f}s; Average processing time: {1:.2f}s; Estimated end of processing in {2:.2f}s)".format(time_taken, average_time, estimated_end))

            pddfall = self.plot_barchart(fusion_frequency_all, "Fusion gene recurrency in single samples", 1)
            pddfall.to_csv(data_out_path_freq_all, sep="\t", index=False)
            #colnames_filtering = ["all fusions", "in/out/neo frame", "both on exon boundary", "with neo peptide", "not black-listed", "by 2-tools", "all together"]
            mmm1 = self.plot_boxswarm(filtering_data_1, colnames_filtering, ["Sample", "Filter", "Data"], "Filter counts of distinct fusion genes - Rep1's")
            fltr_out.write("Rep1's MinMaxMed\tNA\tNA\t" + "\t".join(map(str, mmm1)) + "\n")

            if icam_run:
                pddf12 = self.plot_barchart(fusion_frequency_12, "Fusion gene recurrency in patient samples (replicates merged)", 0)
                pddf12.to_csv(data_out_path_freq_12, sep="\t", index=False)
                mmm2 = self.plot_boxswarm(filtering_data_2, colnames_filtering, ["Sample", "Filter", "Data"], "Filter counts of distinct fusion genes - Rep2's")
                mmm12 = self.plot_boxswarm(filtering_data_12, colnames_filtering, ["Sample", "Filter", "Data"], "Filter counts of distinct fusion genes - Merged")
                fltr_out.write("Rep2's MinMaxMed\tNA\tNA\t" + "\t".join(map(str, mmm2)) + "\n")
                fltr_out.write("Merged MinMaxMed\tNA\tNA\t" + "\t".join(map(str, mmm12)) + "\n")

        with PdfPages(figure_out_path) as pdf:
            for plot_fig in self.figure_list:
                pdf.savefig(plot_fig)

    def plot_barchart(self, fusion_frequency_dict, chart_title, label_cutoff):
        """bkla"""
        freq_df = pd.DataFrame.from_dict(fusion_frequency_dict, orient="index", columns=["Frequency"])
        freq_df["Fusion_Gene"] = freq_df.index
        try:
            plt.figure(len(self.figure_list) + 1)
    #        ax = plt.axes()
            plt.title(chart_title)
            ax = sns.barplot(x="Fusion_Gene", y="Frequency", data=freq_df)
    #        new_x_ticks = []
    #        for i in freq_df.index:
    #            if freq_df.loc[i, "Frequency"] == label_cutoff:
    #                new_x_ticks.append("")
    #            else:
    #                new_x_ticks.append(freq_df.loc[i, "Fusion_Gene"])
            new_x_ticks = []
            for i in freq_df.index:
                if freq_df.loc[i, "Frequency"] == label_cutoff:
                    new_x_ticks.append("")
                else:
                    new_x_ticks.append(freq_df.loc[i, "Fusion_Gene"])
    #        for tick in ax.get_xticklabels():
    #            tick.set_rotation(90)
            ax.set_xticklabels(new_x_ticks, rotation=90)
            plt.tight_layout()
            self.figure_list.append(ax.figure)
        except ValueError:
            print("No fusion genes predicted, i.e. nothing to plot...")
        return freq_df

    def plot_boxswarm(self, data_dict, group_names, column_names, chart_title):
        """Create a pandas dataframe from a dictionary and plot the data table as a heatmap.
        Summarize the data furthermore in a boxplot with swarmplot overlay"""
#        print(data_dict)
#        print(group_names)
        list_of_data_lists = []
        for data_dict_key in data_dict:
            for i, list_record in enumerate(data_dict[data_dict_key]):
                list_of_data_lists.append([data_dict_key, group_names[i], list_record])
        # plot data as box plot (plus swarmplot overlay)
        pd_df = pd.DataFrame(list_of_data_lists, columns=column_names)
        try:
            plt.figure(len(self.figure_list) + 1)
            ax = plt.axes()
            plt.title(chart_title)
            sns.boxplot(data=pd_df, ax=ax, x=column_names[1], y=column_names[2], hue=column_names[1], dodge=False)
            sns.swarmplot(data=pd_df, ax=ax, x=column_names[1], y=column_names[2], color=(0.8, 0.8, 0.8), alpha=0.5)
            ax.legend_.remove()
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)
            plt.tight_layout()
            self.figure_list.append(ax.figure)
        except ValueError:
            print("No fusion genes predicted, i.e. nothing to plot...")

        # get min/max/median for all categories and return them in a list to print
        mmm_list = []
        for group in group_names:
            min_val = pd_df.loc[pd_df[column_names[1]] == group][column_names[2]].min()
            max_val = pd_df.loc[pd_df[column_names[1]] == group][column_names[2]].max()
            med_val = pd_df.loc[pd_df[column_names[1]] == group][column_names[2]].median()
            out_str = str(min_val) + "-" + str(max_val) + "(" + str(med_val) + ")"
            mmm_list.append(out_str)
        return mmm_list

    @staticmethod
    def add_to_fus_dict(input_set, fusion_dict):
        """bla"""
        for fus_gene in input_set:
            if fus_gene not in fusion_dict:
                fusion_dict[fus_gene] = 0
            fusion_dict[fus_gene] += 1
        return fusion_dict

def main():
    """Command line argument parsing and app start"""
    parser = argparse.ArgumentParser(description='Post processing of an easyfuse run - currently, collecting runtimes only :)')

    parser.add_argument('-i', '--input', dest='input', help='Specify the easyfuse root dir of the run you want to process.', required=True)
    parser.add_argument('-c', '--config', dest='config', help='Specify config file.', required=True)
    parser.add_argument('--model_predictions', default=False, action='store_true', help='Score predictions based on pre-train model')
    parser.add_argument('--icam_run', dest='icam_run', help=argparse.SUPPRESS, default=False, action='store_true')
    args = parser.parse_args()
    stats = IcamSummary(args.input, args.config)
    stats.run(args.icam_run, args.model_predictions)

if __name__ == '__main__':
    main()
