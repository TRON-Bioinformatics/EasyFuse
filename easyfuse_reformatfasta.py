"""
Parses an ensembl/gencode derived whole transcriptome multi fasta file
During reading of the fasta file, the header is modified to be in the following format
>"gene_symbol"_"ensembl_gene_id"_"ensembl_transcript_id"
This method needs to be called only once because the modified fasta file will be written
to "~.toNameHeader"
The output of this method is required by the easyfuse_denovoassembly.py script

@author: BNT (URLA)
@version: 20181126
"""

from __future__ import print_function
from argparse import ArgumentParser
import sys

class FastaReformater(object):
    """Select alignments belonging to putative fusions from an s/bam file"""
    def __init__(self, in_fasta):
        self.in_fasta = in_fasta
        self.out_fasta = in_fasta + ".toNameHeader"

    def reformat(self):
        """Parse ensembl transcriptome fasta file and modify header"""
        with open(self.in_fasta, "r") as fasin, open(self.out_fasta, "w") as fasout:
            for line in fasin:
                if line.startswith(">"):
                    try:
                        line_splitter = line.split()
                        transcript_ens_id = line_splitter[0].lstrip(">")
                        gene_ens_id = ""
                        gene_symbol = ""
                        for field in line_splitter:
                            if field.startswith("gene_symbol:"):
                                gene_symbol = field.split(":")[1]
                            elif field.startswith("gene:"):
                                gene_ens_id = field.split(":")[1]
                        fasout.write(">{0}_{1}_{2}\n".format(
                            gene_symbol,
                            gene_ens_id,
                            transcript_ens_id))
                    except IndexError as inderr:
                        print(inderr)
                        print(line)
                        print(field)
                        sys.exit(1)
                else:
                    fasout.write(line)

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input', dest='input', help='Ensembl fasta file', required=True)
    args = parser.parse_args()

    ass = FastaReformater(args.input)
    ass.reformat()

if __name__ == '__main__':
    main()
