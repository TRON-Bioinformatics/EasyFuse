#!/usr/bin/env python3

"""
Get predicted gene fusions from tools
Fusions are parsed into "id / fusion genes / breakpoints / supporting reads / prediction tool"
and output based on their set recurrency value. The class requires the path' of data and
output folders and writes the "Detected_Fusions.csv"

@author: Tron (PASO)
@version: 20241126
"""
from argparse import ArgumentParser
import csv
import logging

from fusionparsing.src.file_headers import OUTPUT_HEADER

FORMAT = '%(asctime)s - %(message)s'
logging.basicConfig(format=FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)


class FusionParser:
    """
    Get and parse results from previously run programs 
    (fusion prediction, hla typing, expression estimation)
    """

    # Initialization of parameters
    def __init__(self, output_path, sample_id):
        """Parameter initiation and work folder creation."""
        self.output_path = output_path
        self.sample_id = sample_id


    def run(self, tool_outputs=None):
        """Parse output of fusion tools and save them into output file"""

        logger.debug("Generating Detected Fusions table")

        count = 0
        with open(self.output_path, "w", encoding="utf8") as all_fusions_file:
            writer = csv.DictWriter(all_fusions_file, fieldnames=OUTPUT_HEADER, delimiter=";")
            writer.writeheader()
            # parses all available fusion tools
            for toolname, tool_output in tool_outputs:
                with open(tool_output, "r", encoding="utf8") as infile:
                    reader = csv.DictReader(infile, delimiter=";")
                    for row in reader:
                        row["Sample"] = self.sample_id
                        row["Tool"] = toolname
                        writer.writerow(row)
                        count += 1

        logger.info("Wrote %s detected fusion genes to %s.", count, self.output_path)


def main():
    """stub"""
    parser = ArgumentParser(description="Parses fusions from specific tools")
    parser.add_argument(
        "--tool",
        dest="input_tools",
        action="append",
        nargs="+",
        metavar=("tool"),
        help="Input file with results from fusion tool (e.g. --tool starfusion results.csv)",
    )

    parser.add_argument(
        "-o",
        "--output_path",
        dest="output_path",
        help="Specify the path to the output file",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--sample",
        dest="sample",
        help="Specify the sample to process.",
        required=True,
    )
    args = parser.parse_args()

    # checks that at least one fusion prediction tool is provided
    assert (args.input_tools), "No fusion prediction tool results provided."

    fusion_parser = FusionParser(
        output_path=args.output_path,
        sample_id=args.sample
    )
    fusion_parser.run(args.input_tools)


if __name__ == "__main__":
    main()
