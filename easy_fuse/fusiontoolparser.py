#!/usr/bin/env python3

"""
Get predicted gene fusions from tools
Fusions are parsed into "id / fusion genes / breakpoints / supporting reads / prediction tool"
and output based on their set recurrency value. The class requires the path' of data and
output folders and writes the "Detected_Fusions.csv"

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""

from argparse import ArgumentParser

import logzero
from logzero import logger

from fusiontoolparser_helper import *


# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class FusionParser(object):
    """Get and parse results from previously run programs (fusion prediction, hla typing, expression estimation)"""

    # Initialization of parameters
    def __init__(self, input_path, output_path, sample_id, tool_num_cutoff, fusiontool_list, sample_log):
        """Parameter initiation and work folder creation."""
        self.input_path = input_path
        self.output_path = output_path
        self.sample_id = sample_id
        self.tool_num_cutoff = int(tool_num_cutoff)
        self.tools = fusiontool_list.split(",")
        logzero.logfile(sample_log)

    def concatenate_fusion_results(self, input_path, output_path):
        """Parse output of fusion tools and save them into output file"""

        count = 0

        logger.debug("Generating Detected Fusions table")
        detected_fusions_file = os.path.join(output_path, "Detected_Fusions.csv")

        with open(detected_fusions_file, "w") as fus_file:
            # write header
            fus_file.write("BPID;Fusion_Gene;Breakpoint1;Breakpoint2;Junction_Reads;Spanning_Reads;Sample;Tool\n")

            for tool in self.tools:
                fusion_map = self.get_tool_results(input_path, output_path, tool)
                for bpid in fusion_map:
                    count += 1
                    fus_file.write("{};{};{};{}\n".format(bpid, ";".join(fusion_map[bpid]), self.sample_id, tool))
        logger.info("Wrote {0} detected fusion genes to {1}.".format(count, detected_fusions_file))

    def get_tool_results(self, input_path, output_path, tool):
        """Return results as dict for a fusion tool"""

        self.logger.info("Parsing results for {}".format(tool))

        fusion_map = parse_results(input_path, tool)

        tool_res_file = os.path.join(output_path, tool.lower() + ".csv")
        with open(tool_res_file, "w") as tool_outf:
            tool_outf.write("bpid;fusion_gene;breakpoint1;breakpoint2;junc_reads;span_reads\n")
            for key in fusion_map:
                tool_outf.write("{};{}\n".format(key, ";".join(fusion_map[key])))
        return fusion_map

    def run(self):
        """This method will parse the results for all individual tools and saves them into a concatenated file"""

        self.concatenate_fusion_results(self.input_path, self.output_path)


def main():
    """Main"""
    parser = ArgumentParser(description='Extracts information on fusion genes')
    parser.add_argument('-i', '--input_path', dest='input_path', help='Specify the folder of the sample.', required=True)
    parser.add_argument('-o', '--output_path', dest='output_path', help='Specify the fusion output folder of the sample.', required=True)
    parser.add_argument('-s', '--sample', dest='sample', help='Specify the sample to process.', required=True)
    parser.add_argument('-t', '--toolcutoff', dest='toolcutoff', help='The minimal of tools that must support a fusion in order to be considered further.', required=True)
    parser.add_argument('-f', '--fusionlist', dest='fusionlist', help='A list of fusion prediction tools that where run on the sample.', required=True)
    parser.add_argument('-l', '--logger', dest='logger', help='Logging of processing steps', default="")
    args = parser.parse_args()

    res_parse_1 = FusionParser(args.input_path, args.output_path, args.sample, args.toolcutoff, args.fusionlist, args.logger)
    res_parse_1.run()
