#!/usr/bin/env python

"""
First truncates the header of a multi fasta file to the first whitespace
Second splits the modified file into separate fasta files
The separated chromosomes are at least mandatory for mapsplice2 and fetchdata
Truncating of the headers is also mandatory for mapsplice2 and in general
a good idea to avoid any kind of parsing issues

@author: BNT (URLA)
@version: 20180131
"""

from __future__ import print_function

import os.path
from argparse import ArgumentParser


class GenomeRefParser(object):
    """Initialize class variables"""

    def __init__(self, in_fasta, out_dir):
        self.in_fasta = in_fasta
        self.out_dir = out_dir

    def short_header(self):
        """Shorten the header in a multi fasta file to the first whitespace"""
        out_fasta = "{}_SH.fa".format(self.in_fasta.rsplit(".", 1)[0])
        with open(self.in_fasta, "r") as fasin, open(out_fasta, "w") as fasout:
            for line in fasin:
                if line.startswith(">"):
                    print("Changing header from \"{0}\" to \"{1}\"".format(
                        line.rstrip("\n"),
                        line.split()[0]))
                    fasout.write("{}\n".format(line.split()[0]))
                else:
                    fasout.write(line)
        self.in_fasta = out_fasta

    def split(self):
        """Performing splitting of multi fasta file"""
        fasout = None
        with open(self.in_fasta, "r") as fasin:
            for line in fasin:
                if line.startswith(">"):
                    if fasout:
                        fasout.close()
                    out_file_name = os.path.join(self.out_dir, "{}.fa".format(
                        line.split()[0][1:]))
                    print("Writing fasta record \"{0}\" to file \"{1}\"".format(
                        line.rstrip("\n"),
                        out_file_name))
                    fasout = open(out_file_name, "w")
                fasout.write(line)
            if fasout:
                fasout.close()


def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input', dest='input', help='multi fasta file', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Output directory', required=True)
    args = parser.parse_args()

    grp = GenomeRefParser(args.input, args.output)
    grp.short_header()
    grp.split()


if __name__ == '__main__':
    main()
