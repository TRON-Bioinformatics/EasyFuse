#!/usr/bin/env python3

"""
@author: BNT (URLA), TRON (PASO)
@version: 20190118
"""
import os
import os.path
import time
from typing import List

from logzero import logger

import easy_fuse.misc.io_methods as IOMethods
from easy_fuse.join_data import DataJoining
from easy_fuse.misc.config import EasyFuseConfiguration


class FusionSummary(object):
    """Collect stats of the run and write them to file"""

    def __init__(self, input_path, config: EasyFuseConfiguration, samples: List[str]):
        self.input_path = input_path
        self.cfg = config
        self.samples = samples

    def run(self, model_predictions):
        """Execute individual methods"""
        fusion_data_summary_path = os.path.join(self.input_path, "FusionSummary")
        IOMethods.create_folder(fusion_data_summary_path)

        fusion_frequency_all = {}
        filtering_data_1 = {}
        # 0:type; 1:boundary; 2:frame; 3:pepseq; 4:counts; 5:blacklist; 6:prediction;
        # 7:exonError; 8:unfiltered; 9:allFilter; 10:allButPredFilter

        average_time = 0
        count_processed = 0

        for sample in self.samples:
            count_processed += 1
            print("Processing sample {0} (dataset {1}/{2})".format(sample, count_processed, len(self.samples)))
            start_time = time.time()
            fusion_data_summary = DataJoining(self.input_path, sample, "",
                                              os.path.join(fusion_data_summary_path, sample), model_predictions,
                                              self.cfg).run()
            fusion_frequency_all = self.add_to_fus_dict(fusion_data_summary[1], fusion_frequency_all)
            filtering_data_1[sample] = fusion_data_summary[0]

            time_taken = time.time() - start_time
            average_time = (average_time * (count_processed - 1) + time_taken) / count_processed
            logger.info(
                "Done. Processing time: {0:.2f}s; Average processing time: {1:.2f}s)\n".format(
                    time_taken, average_time))

    @staticmethod
    def add_to_fus_dict(input_set, fusion_dict):
        """bla"""
        for fus_gene in input_set:
            if fus_gene not in fusion_dict:
                fusion_dict[fus_gene] = 0
            fusion_dict[fus_gene] += 1
        return fusion_dict


def add_summarize_data_args(parser):
    parser.add_argument('-i', '--input', dest='input',
                        help='Specify the easyfuse root dir of the run you want to process.', required=True)
    parser.add_argument('-c', '--config-file', dest='config_file', required=True,
                        help='Specify alternative config file to use for your analysis')

    parser.add_argument('--model_predictions', default=False, action='store_true',
                        help='Score predictions based on pre-train model')
    parser.add_argument('--samples', default=False, nargs='+',
                        help='List of samples to analyze')
    parser.set_defaults(func=summarize_data_command)


def summarize_data_command(args):
    config = EasyFuseConfiguration(args.config_file)
    summary = FusionSummary(args.input, config, args.samples)
    summary.run(args.model_predictions)
