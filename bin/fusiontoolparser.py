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
import os

FORMAT = '%(asctime)s - %(message)s'
logging.basicConfig(format=FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)


from fusiontoolparser.file_headers import TOOL_HEADER, OUTPUT_HEADER
from fusiontoolparser.fusiontoolparser_helper import parse_fusioncatcher_results
from fusiontoolparser.fusiontoolparser_helper import parse_arriba_results
from fusiontoolparser.fusiontoolparser_helper import parse_mapsplice_results
from fusiontoolparser.fusiontoolparser_helper import parse_starfusion_results
from fusiontoolparser.fusiontoolparser_helper import parse_soapfuse_results
from fusiontoolparser.fusiontoolparser_helper import parse_infusion_results


class FusionParser:
    """
    Get and parse results from previously run programs 
    (fusion prediction, hla typing, expression estimation)
    """

    # Initialization of parameters
    def __init__(
            self, input_fusioncatcher, input_fusioncatcher2, input_infusion, input_mapsplice,
            input_starfusion, input_soapfuse, input_arriba,
            output_path, sample_id):
        """Parameter initiation and work folder creation."""
        self.input_fusioncatcher = input_fusioncatcher
        self.input_fusioncatcher2 = input_fusioncatcher2
        self.input_infusion = input_infusion
        self.input_mapsplice = input_mapsplice
        self.input_starfusion = input_starfusion
        self.input_soapfuse = input_soapfuse
        self.input_arriba = input_arriba
        self.output_path = output_path
        self.sample_id = sample_id


    def run(self):
        """Parse output of fusion tools and save them into output file"""

        count = 0

        logger.debug("Generating Detected Fusions table")
        detected_fusions_file = os.path.join(self.output_path, "Detected_Fusions.csv")

        with open(detected_fusions_file, "w", encoding="utf8") as all_fusions_file:
            writer = csv.writer(all_fusions_file, fieldnames=OUTPUT_HEADER, delimiter=";")
            writer.writeheader()
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
            if self.input_arriba and os.path.exists(self.input_arriba):
                fusions = parse_arriba_results(self.input_arriba)
                count += self._write_fusions(
                    all_fusions_file=all_fusions_file, fusions_dict=fusions, tool="arriba")

        logger.info("Wrote %s detected fusion genes to %s.", count, detected_fusions_file)

    def _write_fusions(self, fusions_dict: dict, all_fusions_file, tool):
        count = 0

        # write fusions for this tool to a file
        tool_res_file = os.path.join(self.output_path, tool.lower() + ".csv")
        with open(tool_res_file, "w", encoding="utf8") as tool_outf:
            writer = csv.writer(tool_outf, fieldnames=TOOL_HEADER, delimiter=";")
            writer.writeheader()
            for key in fusions_dict:
                writer.writerow(fusions_dict[key])

        # add fusions from this tool to all fusions file
        for bpid in fusions_dict:
            count += 1
            all_fusions_file.write(
                "{};{};{};{}\n".format(
                    bpid, ";".join(fusions_dict[bpid]), self.sample_id, tool
                )
            )
        return count


def main():
    """stub"""
    parser = ArgumentParser(description="Parses fusions from specific tools")
    parser.add_argument(
        "--input_fusioncatcher",
        dest="input_fusioncatcher",
        help="Input file with results from fusion catcher (ie: summary_candidate_fusions.txt).",
    )
    parser.add_argument(
        "--input_fusioncatcher2",
        dest="input_fusioncatcher2",
        help="Input file with results from fusion catcher (ie: final-list_candidate-fusion-genes.txt).",
    )
    parser.add_argument(
        "--input_starfusion",
        dest="input_starfusion",
        help="Input file with results from star fusion (ie: star-fusion.fusion_predictions.tsv).",
    )
    parser.add_argument(
        "--input_infusion",
        dest="input_infusion",
        help="Input file with results from infusion (ie: final_fusion_genes).",
    )
    parser.add_argument(
        "--input_mapsplice",
        dest="input_mapsplice",
        help="Input file with results from mapsplice (ie: fusions_well_annotated.txt).",
    )
    parser.add_argument(
        "--input_soapfuse",
        dest="input_soapfuse",
        help="Input file with results from soapfuse "
             "(ie: final_fusion_genes/${SAMPLE}/${SAMPLE}.final.Fusion.specific.for.genes).",
    )
    parser.add_argument(
        "--input_arriba",
        dest="input_arriba",
        help="Input file with results from arriba (ie: fusions.tsv).",
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
    args = parser.parse_args()

    # checks that at least one fusion prediction tool is provided
    assert (args.input_fusioncatcher or
            args.input_starfusion or
            args.input_infusion or
            args.input_mapsplice or
            args.input_soapfuse or
            args.input_arriba), "No fusion prediction tool results provided."

    fusion_parser = FusionParser(
        input_fusioncatcher=args.input_fusioncatcher,
        input_fusioncatcher2=args.input_fusioncatcher2,
        input_starfusion=args.input_starfusion,
        input_infusion=args.input_infusion,
        input_mapsplice=args.input_mapsplice,
        input_soapfuse=args.input_soapfuse,
        input_arriba=args.input_arriba,
        output_path=args.output_path,
        sample_id=args.sample
    )
    fusion_parser.run()


if __name__ == "__main__":
    main()
