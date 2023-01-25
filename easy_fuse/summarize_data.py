#!/usr/bin/env python3

"""
@author: BNT (URLA), TRON (PASO)
@version: 20190118
"""
import os
import os.path
import time
from argparse import ArgumentParser

import easy_fuse.misc.io_methods as IOMethods
from easy_fuse.join_data import DataJoining
from easy_fuse.misc.config import EasyFuseConfiguration
from easy_fuse.misc.samples import SamplesDB


class FusionSummary(object):
    """Collect stats of the run and write them to file"""

    def __init__(self, input_path, config: EasyFuseConfiguration):
        self.input_path = input_path
        self.samples = SamplesDB(os.path.join(input_path, "samples.db"))
        self.cfg = config

    def run(self, model_predictions):
        """Execute individual methods"""
        fusion_tools = self.cfg["general"]["fusiontools"].split(",")
        fusion_data_summary_path = os.path.join(self.input_path, "FusionSummary")
        IOMethods.create_folder(fusion_data_summary_path)

        fusion_frequency_all = {}
        filtering_data_1 = {}
        # 0:type; 1:boundary; 2:frame; 3:pepseq; 4:counts; 5:blacklist; 6:prediction;
        # 7:exonError; 8:unfiltered; 9:allFilter; 10:allButPredFilter
        # get sorted list of sample ids
        sid_list = sorted(self.samples.get_sample_ids())
        #        sid_list = sid_list[:10]

        count_valid_sample = 0
        # check how many samples are available for processing and create {sample_id: sample_date} dict
        sample_date_dict = {}
        sample_toolCnt_dict = {}
        for i, sample in enumerate(sid_list, 1):
            #            try:
            #                sample_date_dict[sample] = self.samples.get_fastq_files(sample)[0].split("/")[8].split("_")[0]
            #            except IndexError:
            #                print("Warning: Could not retrieve sample date for {}".format(sample))
            sample_date_dict[sample] = "NA"
            # count the numbers of fusion tools that have been run on this sample
            sample_toolCnt_dict[sample] = len(
                [tool for tool in fusion_tools if tool in self.samples.get_tool_list_from_state(sample)])
            if "fetchdata" in self.samples.get_tool_list_from_state(sample):
                count_valid_sample += 1
        print(
            "Found {0} (partially) processed samples in {1}. Data will be collected from {2} samples for which fetchdata has been run.".format(
                i, self.input_path, count_valid_sample))

        average_time = 0
        count_processed = 0

        for sample in sid_list:

            if "requantifyBest" in self.samples.get_tool_list_from_state(sample):
                count_processed += 1
                print("Processing sample {0} (dataset {1}/{2})".format(sample, count_processed, len(sid_list)))
                start_time = time.time()
                fusion_data_summary = DataJoining(self.input_path, sample, "",
                                                  os.path.join(fusion_data_summary_path, sample), model_predictions,
                                                  self.cfg).run()
                fusion_frequency_all = self.add_to_fus_dict(fusion_data_summary[1], fusion_frequency_all)
                filtering_data_1[sample] = fusion_data_summary[0]

                time_taken = time.time() - start_time
                average_time = (average_time * (count_processed - 1) + time_taken) / count_processed
                estimated_end = average_time * (count_valid_sample - count_processed)
                print(
                    "done. Processing time: {0:.2f}s; Average processing time: {1:.2f}s; Estimated end of processing in {2:.2f}s)\n".format(
                        time_taken, average_time, estimated_end))

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
    parser = ArgumentParser(description='Post processing of an easyfuse run - currently, collecting runtimes only :)')

    parser.add_argument('-i', '--input', dest='input',
                        help='Specify the easyfuse root dir of the run you want to process.', required=True)
    parser.add_argument('-c', '--config-file', dest='config_file', required=True,
                        help='Specify alternative config file to use for your analysis')

    parser.add_argument('--model_predictions', default=False, action='store_true',
                        help='Score predictions based on pre-train model')
    args = parser.parse_args()
    config = EasyFuseConfiguration(args.config_file)
    stats = FusionSummary(args.input, config)
    stats.run(args.model_predictions)
