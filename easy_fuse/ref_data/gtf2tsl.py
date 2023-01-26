#!/usr/bin/env python

"""
Get Transcript Support Level (TSL) info from an ensembl GTF file

@author: BNT (URLA)
@version: 20190702
"""

import argparse
import sys


class Gtf2Tsl(object):
    """stub"""

    def __init__(self, gtf_in, tsl_out):
        """Parameter initialization"""
        self.gtf_in = gtf_in
        self.tsl_out = tsl_out

    def get_data(self):
        """ stub """
        feature_list = []
        tsl_dict = {}
        count_transcripts = 0
        with open(self.gtf_in, "r") as in_file:
            for line in in_file:
                if line.startswith("#"):  # skip header
                    continue
                _, _, feature, _, _, _, _, _, attribute = line.split("\t")
                if not feature in feature_list:
                    feature_list.append(feature)
                if feature == "transcript":
                    count_transcripts += 1
                    field_split = attribute.split(";")
                    transcript_id = "NA"
                    trans_biotype = "NA"
                    gene_biotype = "NA"
                    tsl = "NA"
                    for field in field_split:
                        transcript_id = self.get_field_value(transcript_id, field, "transcript_id")
                        trans_biotype = self.get_field_value(trans_biotype, field, "transcript_biotype")
                        gene_biotype = self.get_field_value(gene_biotype, field, "gene_biotype")
                        tsl = self.get_field_value(tsl, field, "transcript_support_level")

                    print("{}_{}_{}_{}".format(transcript_id, trans_biotype, gene_biotype, tsl))
                    # check again, that what we have is what we expect
                    if not (transcript_id[0:3] == "ENS" or transcript_id[0:4] == "MGP_"):
                        print(
                            "This does not look like a valid ensembl id: {}. Is this a valid ensembl gtf line: {}?".format(
                                transcript_id, line))
                        sys.exit(99)
                    if not transcript_id in tsl_dict:
                        tsl_dict[transcript_id] = [trans_biotype, gene_biotype, tsl]
                    if count_transcripts == 20:
                        break

        # write data to file
        with open(self.tsl_out, "w") as out_file:
            out_file.write("{}\n".format("\t".join(["transcript_id", "trans_biotype", "gene_biotype", "tsl"])))
            for transcript_id in tsl_dict:
                out_file.write("{}\t{}\n".format(transcript_id, "\t".join(tsl_dict[transcript_id])))
        # list of available feature strings
        print(feature_list)

    def get_field_value(self, old_value, field_string, lookup_string):
        """ Helper method to get the value from a key-value pair in the form of: ( this the value "and this the key") """
        if not old_value == "NA":  # check if value has been found previously
            return old_value
        new_value = "NA"
        field_string = field_string.strip()  # remove leading/trailing whitespaces
        if field_string.startswith(lookup_string):
            new_value = field_string.split("\"")[1]  # get the value inside quotation marks
        return new_value


def main():
    """Parse command line arguments and start script"""
    parser = argparse.ArgumentParser(description="Get Transcript Support Level (TSL) info from an ensembl GTF file")
    parser.add_argument('--gtf', dest='gtf', help='Gtf input file', required=True)
    parser.add_argument('--tsl', dest='tsl', help='Tsl output file', required=True)
    args = parser.parse_args()

    Gtf2Tsl(args.gtf, args.tsl).get_data()


if __name__ == '__main__':
    main()
