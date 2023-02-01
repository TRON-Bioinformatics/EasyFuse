#!/usr/bin/env python3

"""
Get predicted gene fusions from tools
Fusions are parsed into "id / fusion genes / breakpoints / supporting reads / prediction tool"
and output based on their set recurrency value. The class requires the path' of data and
output folders and writes the "Detected_Fusions.csv"

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""
import os
from argparse import ArgumentParser
from typing import Dict

import logzero
from logzero import logger

from easy_fuse.fusiontoolparser_helper import *


# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class FusionParser(object):
    """Get and parse results from previously run programs (fusion prediction, hla typing, expression estimation)"""

    # Initialization of parameters
    def __init__(
            self, input_fusioncatcher, input_fusioncatcher2, input_infusion, input_mapsplice,
            input_starfusion, input_soapfuse,
            output_path, sample_id, logfile):
        """Parameter initiation and work folder creation."""
        self.input_fusioncatcher = input_fusioncatcher
        self.input_fusioncatcher2 = input_fusioncatcher2
        self.input_infusion = input_infusion
        self.input_mapsplice = input_mapsplice
        self.input_starfusion = input_starfusion
        self.input_soapfuse = input_soapfuse
        self.output_path = output_path
        self.sample_id = sample_id
        logzero.logfile(logfile)

    def run(self):
        """Parse output of fusion tools and save them into output file"""

        count = 0

        logger.debug("Generating Detected Fusions table")
        detected_fusions_file = os.path.join(self.output_path, "Detected_Fusions.csv")

        with open(detected_fusions_file, "w") as all_fusions_file:
            # write header
            all_fusions_file.write(
                "BPID;Fusion_Gene;Breakpoint1;Breakpoint2;Junction_Reads;Spanning_Reads;Sample;Tool\n"
            )
            # parses all available fusion tools
            if self.input_fusioncatcher and os.path.exists(self.input_fusioncatcher):
                fusions = parse_fusioncatcher_results(self.input_fusioncatcher, self.input_fusioncatcher2)
                count += self._write_fusions(
                    all_fusions_file=all_fusions_file, fusions_dict=fusions, tool="fusioncatcher")
            if self.input_infusion and os.path.exists(self.input_infusion):
                fusions = parse_infusion_results(self.input_infusion)
                count += self._write_fusions(
                    all_fusions_file=all_fusions_file, fusions_dict=fusions, tool="infusion")
            if self.input_mapsplice and os.path.exists(self.input_mapsplice):
                fusions = parse_mapsplice_results(self.input_mapsplice)
                count += self._write_fusions(
                    all_fusions_file=all_fusions_file, fusions_dict=fusions, tool="mapsplice")
            if self.input_starfusion and os.path.exists(self.input_starfusion):
                fusions = parse_starfusion_results(self.input_starfusion)
                count += self._write_fusions(
                    all_fusions_file=all_fusions_file, fusions_dict=fusions, tool="starfusion")
            if self.input_soapfuse and os.path.exists(self.input_soapfuse):
                fusions = parse_soapfuse_results(self.input_soapfuse)
                count += self._write_fusions(
                    all_fusions_file=all_fusions_file, fusions_dict=fusions, tool="soapfuse")

        logger.info(
            "Wrote {0} detected fusion genes to {1}.".format(
                count, detected_fusions_file
            )
        )

    def _write_fusions(self, fusions_dict: Dict, all_fusions_file, tool):
        count = 0

        # write fusions for this tool to a file
        tool_res_file = os.path.join(self.output_path, tool.lower() + ".csv")
        with open(tool_res_file, "w") as tool_outf:
            tool_outf.write(
                "bpid;fusion_gene;breakpoint1;breakpoint2;junc_reads;span_reads\n"
            )
            for key in fusions_dict:
                tool_outf.write("{};{}\n".format(key, ";".join(fusions_dict[key])))

        # add fusions from this tool to all fusions file
        for bpid in fusions_dict:
            count += 1
            all_fusions_file.write(
                "{};{};{};{}\n".format(
                    bpid, ";".join(fusions_dict[bpid]), self.sample_id, tool
                )
            )
        return count


def add_fusion_parse_args(parser):

    parser.add_argument(
        "--input-fusioncatcher",
        dest="input_fusioncatcher",
        help="Input file with results from fusion catcher (ie: summary_candidate_fusions.txt).",
    )
    parser.add_argument(
        "--input-fusioncatcher2",
        dest="input_fusioncatcher2",
        help="Input file with results from fusion catcher (ie: final-list_candidate-fusion-genes.txt).",
    )
    parser.add_argument(
        "--input-starfusion",
        dest="input_starfusion",
        help="Input file with results from star fusion (ie: star-fusion.fusion_predictions.tsv).",
    )
    parser.add_argument(
        "--input-infusion",
        dest="input_infusion",
        help="Input file with results from infusion (ie: final_fusion_genes).",
    )
    parser.add_argument(
        "--input-mapsplice",
        dest="input_mapsplice",
        help="Input file with results from mapsplice (ie: fusions_well_annotated.txt).",
    )
    parser.add_argument(
        "--input-soapfuse",
        dest="input_soapfuse",
        help="Input file with results from soapfuse "
             "(ie: final_fusion_genes/${SAMPLE}/${SAMPLE}.final.Fusion.specific.for.genes).",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_path",
        help="Specify the fusion output folder of the sample.",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--sample",
        dest="sample",
        help="Specify the sample to process.",
        required=True,
    )
    parser.add_argument(
        "-l", "--logger", dest="logger", help="Logging of processing steps", default=""
    )
    parser.set_defaults(func=fusion_parser_command)


def fusion_parser_command(args):

    # checks that at least one fusion prediction tool is provided
    assert args.input_fusioncatcher or args.input_starfusion or args.input_infusion or args.input_mapsplice or \
           args.input_soapfuse, "No fusion prediction tool results provided."

    fusion_parser = FusionParser(
        input_fusioncatcher=args.input_fusioncatcher,
        input_fusioncatcher2=args.input_fusioncatcher2,
        input_starfusion=args.input_starfusion,
        input_infusion=args.input_infusion,
        input_mapsplice=args.input_mapsplice,
        input_soapfuse=args.input_soapfuse,
        output_path=args.output_path,
        sample_id=args.sample,
        logfile=args.logger
    )
    fusion_parser.run()
