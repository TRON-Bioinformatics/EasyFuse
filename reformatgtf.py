"""
Parses an ensembl/gencode derived gene annotation file (gtf)
During reading of the gtf file, the coordinates of all "genes" are exstracted
and written to a new tdt file in the following format
"gene_symbol"\t"chr"\t"start"\t"stop"
This method needs to be called only once because the modified fasta file will be written
to "~.togenecoords"
The output of this method is required by the easyfuse_denovoassembly.py script

@author: BNT (URLA)
@version: 20181126
"""

from __future__ import print_function
from argparse import ArgumentParser

class GtfReformater(object):
    """Select alignments belonging to putative fusions from an s/bam file"""
    def __init__(self, in_gtf):
        self.in_gtf = in_gtf
        self.out_gtf = in_gtf + ".togenecoords"

    def reformat(self):
        """change gtf"""
        with open(self.in_gtf, "r") as gtfin, open(self.out_gtf, "w") as gtfout:
            for line in gtfin:
                if line.startswith("#"):
                    continue
                line_splitter = line.split("\t")
                if line_splitter[2] == "gene":
                    field_splitter = line_splitter[8].split(";")
                    for field in field_splitter:
                        if "gene_name" in field:
                            gtfout.write(
                                "\t".join([
                                    field.lstrip(" gene_name ").strip("\""),
                                    line_splitter[0],
                                    line_splitter[3],
                                    line_splitter[4]
                                    ])
                                + "\n"
                                )

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input', dest='input', help='Ensembl GTF file', required=True)
    args = parser.parse_args()

    ass = GtfReformater(args.input)
    ass.reformat()

if __name__ == '__main__':
    main()
