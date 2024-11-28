#!/usr/bin/env python3

"""
This module merges the results from the individual fusion tools
into a single output file

@author: Tron (PASO)
@version: 20241127
"""

from argparse import ArgumentParser
import csv

from fusionparsing.src.file_headers import TOOL_HEADER
from fusionparsing.src.parse_arriba import parse_arriba_results
from fusionparsing.src.parse_fusioncatcher import parse_fusioncatcher_results
from fusionparsing.src.parse_infusion import parse_infusion_results
from fusionparsing.src.parse_mapsplice import parse_mapsplice_results
from fusionparsing.src.parse_starfusion import parse_starfusion_results
from fusionparsing.src.parse_soapfuse import parse_soapfuse_results


def parse_tool_results(infile: str, infile2: str, tool: str) -> list:
    """Returns the results of the specified tool

    Args:
        infile (str): Path to first input result file
        infile2 (str): Path to second optional input result file
        tool (str): Tool name

    Returns:
        list: List of fusion results
    """
    result = []
    if tool == "arriba":
        result = parse_arriba_results(infile)
    if tool == "fusioncatcher":
        result = parse_fusioncatcher_results(infile, infile2)
    if tool == "infusion":
        result = parse_infusion_results(infile)
    if tool == "mapsplice":
        result = parse_mapsplice_results(infile)
    if tool == "starfusion":
        result = parse_starfusion_results(infile)
    if tool == "soapfuse":
        result = parse_soapfuse_results(infile)
    return result


def write_results(fusions: list, output_file: str):
    """Write the fusion results to the output file

    Args:
        fusions (list): List of fusion events predicted by the tool
        output_file (str): Path to the output file
    """
    with open(output_file, "w", encoding="utf8") as infile:
        writer = csv.DictWriter(infile, fieldnames=TOOL_HEADER, delimiter=";")
        writer.writeheader()
        for fusion in fusions:
            writer.writerow(fusion)


def run(input_file: str, input_file2: str, output_file: str, tool: str):
    """Runs the fusion tool parser

    Args:
        input_file (str): Path to first input file
        input_file2 (str): Path to second input file, if required
        output_file (str): Path to the output file
        tool (str): Name of the fusion tool
    """
    fusion_map = parse_tool_results(input_file, input_file2, tool)
    write_results(fusion_map, output_file)


def main():
    """main stub"""
    parser = ArgumentParser(description="Get predicted gene fusions from tools")
    parser.add_argument(
        "-i", "--input_file",
        dest="input_file",
        help="Input file",
        required=True
    )
    parser.add_argument(
        "-I", "--input_file2",
        dest="input_file2",
        help="Input file 2",
        default=""
    )
    parser.add_argument(
        "-o", "--output_file",
        dest="output_file",
        help="Output file",
        required=True
    )
    parser.add_argument(
        "-t", "--tool",
        dest="tool",
        help="Tool name",
        required=True
    )

    args = parser.parse_args()
    run(args.input_file, args.input_file2, args.output_file, args.tool)


if __name__ == "__main__":
    main()
