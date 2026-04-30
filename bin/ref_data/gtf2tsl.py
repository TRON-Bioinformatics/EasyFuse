#!/usr/bin/env python

"""
Get Transcript Support Level (TSL) info from an ensembl GTF file

@author: TRON (PASO)
@version: 20241119
"""

from argparse import ArgumentParser
import csv
import logging
import sys

# FORMAT = '%(asctime)s - %(message)s'
# logging.basicConfig(format=FORMAT, level=logging.INFO)
# logger = logging.getLogger(__name__)


ENS_PREFIX = "ENS"
MGP_PREFIX = "MGP_"

TSL_HEADER = [
    "transcript_id",
    "trans_biotype",
    "gene_biotype",
    "tsl"
]


class Gtf2Tsl:
    """Class for extracting TSL info from an ensembl GTF file"""
    def __init__(self, gtf_in, tsl_out):
        """Parameter initialization"""
        self.gtf_in = gtf_in
        self.tsl_out = tsl_out


    def get_field_value(self, old_value, field_string, lookup_string):
        """
        Helper method to get the value from a key-value pair
        in the form of: ( this the value "and this the key")
        """

        if old_value:
            return old_value
        new_value = "NA"
        field_string = field_string.strip()
        if field_string.startswith(lookup_string):
            new_value = field_string.split("\"")[1]
        return new_value


    def get_tsl(self, attributes: str) -> tuple:
        """Get attribute fields from GTF file.

        Args:
            attributes (str): Attributes field from GTF file

        Returns:
            tuple: Fields from attributes field
        """
        field_split = attributes.split(";")
        transcript_id = ""
        trans_biotype = ""
        gene_biotype = ""
        tsl = ""
        for field in field_split:
            transcript_id = self.get_field_value(
                transcript_id, field, "transcript_id"
            )
            trans_biotype = self.get_field_value(
                trans_biotype, field, "transcript_biotype"
            )
            gene_biotype = self.get_field_value(
                gene_biotype, field, "gene_biotype"
            )
            tsl = self.get_field_value(
                tsl, field, "transcript_support_level"
            )
        return (transcript_id, trans_biotype, gene_biotype, tsl)


    def check_prefix(self, transcript_id: str, row: dict):
        """Check if transcript id has the correct prefix"""
        if not (transcript_id[0:3] == ENS_PREFIX or transcript_id[0:4] == MGP_PREFIX):
            # logger.info("This does not look like a valid ensembl id: %s", transcript_id)
            # logger.info("Is this a valid ensembl gtf line: %s?", row)
            sys.exit(99)


    def parse_gtf(self) -> tuple:
        """Parse GTF file and extract TSL info"""
        feature_list = []
        tsl_dict = {}
        count_transcripts = 0
        with open(self.gtf_in, "r", encoding="utf8") as in_file:
            csvreader = csv.DictReader(filter(lambda row: row[0]!='#', in_file), delimiter="\t")
            for row in csvreader:
                if not row["feature"] in feature_list:
                    feature_list.append(row["feature"])

                if row["feature"] == "transcript":
                    count_transcripts += 1
                    (transcript_id, trans_biotype, gene_biotype, tsl) = self.get_tsl(
                        row["attributes"]
                    )
                    info_str = f"{transcript_id}_{trans_biotype}_{gene_biotype}_{tsl}"
                    # logger.info(info_str)
                    # check again, that what we have is what we expect
                    self.check_prefix(transcript_id, row)
                    if not transcript_id in tsl_dict:
                        tsl_dict[transcript_id] = [trans_biotype, gene_biotype, tsl]
                    if count_transcripts == 20:
                        break
        return tsl_dict, feature_list


    def write_outfile(self, tsl_dict: dict):
        """Write TSL data to file."""
        # logger.info("Writing TSL data to %s", self.tsl_out)
        with open(self.tsl_out, "w", encoding="utf8") as out_file:
            header_str = "\t".join(TSL_HEADER)
            out_file.write(f"{header_str}\n")
            for transcript_id, tsl_info in tsl_dict.items():
                tsl_info_str = "\t".join(tsl_info)
                out_file.write(f"{transcript_id}\t{tsl_info_str}\n")


    def run(self):
        """ stub """
        tsl_dict, feature_list = self.parse_gtf()
        self.write_outfile(tsl_dict)
        # list of available feature strings
        # logger.info(feature_list)


def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(
        description="Get Transcript Support Level (TSL) info from an ensembl GTF file"
    )
    parser.add_argument(
        '--gtf',
        dest='gtf',
        help='Gtf input file',
        required=True
    )
    parser.add_argument(
        '--tsl',
        dest='tsl',
        help='Tsl output file',
        required=True
    )
    args = parser.parse_args()

    Gtf2Tsl(args.gtf, args.tsl).run()


if __name__ == '__main__':
    main()
