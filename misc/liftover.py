#!/usr/bin/env python

"""
Changes the coordinates from one genome annotation to another (liftover) with crossmap
Coordinates are changed in DetectedFusions.csv; all intermediate files are preserved.

@author: BNT (URLA)
@version: 20190131
"""

from __future__ import print_function
from argparse import ArgumentParser
import sys
import os.path
import pandas as pd
from shutil import copy2
import queue as Queueing
from config import Config
from logger import Logger

class FusionLiftover(object):
    """Select alignments belonging to putative fusions from an s/bam file"""
    def __init__(self, in_fus_detect, config, logger):
        self.in_fus_detect = in_fus_detect
        self.cfg = config
        self.logger = Logger(logger)
        self.chr_list = [
            "1", "2", "3", "4", "5",
            "6", "7", "8", "9", "10",
            "11", "12", "13", "14", "15",
            "16", "17", "18", "19", "20",
            "21", "22", "X", "Y", "MT"
        ]

        copy2(in_fus_detect, "{}.bak".format(in_fus_detect))

    def liftcoords(self):
        """Parse ensembl transcriptome fasta file and modify header"""
        crossmap_chain = self.cfg.get('liftover', 'crossmap_chain')

        # check whether annotations and references are set, existing and matching
        ref_genome_org = os.path.basename(crossmap_chain).replace(".", "_").split("_")[0]
        ref_genome_dest = os.path.basename(crossmap_chain).replace(".", "_").split("_")[2]
        if not ref_genome_org.lower() == self.cfg.get('general', 'ref_genome_build').lower():
            self.logger.error("Error 99: Mismatch between set genome version and chain file. Please verify that you have the correct chain file was supplied!")
            print("Error 99: Mismatch between set genome version and chain file. Please verify that you have the correct chain file was supplied!")
            sys.exit(99)

        self.logger.info("Performing liftover from {} to {} with crossmap".format(ref_genome_org, ref_genome_dest))
        tmp_org_bed = "{}.liftover_{}.bed".format(self.in_fus_detect[:-4], ref_genome_org)
        tmp_dest_bed = "{}.liftover_{}.bed".format(self.in_fus_detect[:-4], ref_genome_dest)
        tmp_dest_bed_unmap = "{}.liftover_{}.bed.unmap".format(self.in_fus_detect[:-4], ref_genome_dest)

        # read detected fusions and create bed file from bp coords
        # urla - note: I store the input data in a pandas df as the lifted data can be easily integrated, checks be performed and results written
        in_fus_detect_pddf = pd.read_csv(self.in_fus_detect, sep=";")
        with open(tmp_org_bed, "w") as fusout:
            for i in in_fus_detect_pddf.index:
                chrom, pos, strand = in_fus_detect_pddf.loc[i, "Breakpoint1"].split(":")
                fusout.write("{0}\t{1}\t{2}\t{3}\t1\t{4}\n".format(chrom, pos, (int(pos) + 1), strand, in_fus_detect_pddf.loc[i, "FGID"]))
                chrom, pos, _ = in_fus_detect_pddf.loc[i, "Breakpoint2"].split(":")
                fusout.write("{0}\t{1}\t{2}\t{3}\t2\t{4}\n".format(chrom, pos, (int(pos) + 1), strand, in_fus_detect_pddf.loc[i, "FGID"]))

        # run crossmap to perform liftover
        cmd_crossmap = "{0} bed {1} {2} {3}".format(self.cfg.get('commands', 'crossmap_cmd'), crossmap_chain, tmp_org_bed, tmp_dest_bed)
        Queueing.submit("", cmd_crossmap.split(" "), "", "", "", "", "", "", "", "", "none")
        # check whether some coords were unmapped and print which those are (i.e. which fusion will be lost)
        if os.stat(tmp_dest_bed_unmap).st_size == 0:
            self.logger.info("Liftover was successful for all fusion breakpoints! Creating a new DetectedFusions table...")
        else:
            self.logger.debug("Fusions couldn't be lifted completely. The following fusions have to be excluded from downstream analyses.")
            with open(tmp_dest_bed_unmap, "r") as liftout:
                for line in liftout:
                    print(line)

        # convert lifted bed back to a proper detected fusions table
        in_fus_detect_pddf["bp1_lifted"] = ""
        in_fus_detect_pddf["bp2_lifted"] = ""
        # for each line in the bed file, lookup its corresponding record in the detected fusions file,
        # re-construct the breakpoint string
        # and save them in the additional columns
        # urla - note: probably not the fastest way, but the safest :)

        with open(tmp_dest_bed, "r") as liftout:
            for line in liftout:
#                print(line.rstrip())
                line_splitter = line.rstrip("\n").split("\t")
                if line_splitter[0] in self.chr_list:
                    if line_splitter[4] == "1":
                        in_fus_detect_pddf.loc[in_fus_detect_pddf["FGID"] == line_splitter[5], "bp1_lifted"] = ":".join([line_splitter[0], line_splitter[1], line_splitter[3]])
                    elif line_splitter[4] == "2":
                        in_fus_detect_pddf.loc[in_fus_detect_pddf["FGID"] == line_splitter[5], "bp2_lifted"] = ":".join([line_splitter[0], line_splitter[1], line_splitter[3]])

        for i in in_fus_detect_pddf.index:
            # replace the old breakpoints, with the new ones in the FGID
            if in_fus_detect_pddf.loc[i, "bp1_lifted"] and in_fus_detect_pddf.loc[i, "bp2_lifted"]:
#                print(in_fus_detect_pddf.loc[i, "bp1_lifted"])
#                print(in_fus_detect_pddf.loc[i, "bp2_lifted"])
                in_fus_detect_pddf.loc[i, "FGID"] = in_fus_detect_pddf.loc[i, "FGID"].replace(in_fus_detect_pddf.loc[i, "Breakpoint1"], in_fus_detect_pddf.loc[i, "bp1_lifted"])
                in_fus_detect_pddf.loc[i, "FGID"] = in_fus_detect_pddf.loc[i, "FGID"].replace(in_fus_detect_pddf.loc[i, "Breakpoint2"], in_fus_detect_pddf.loc[i, "bp2_lifted"])
                # then overwrite the old breakpoints
                in_fus_detect_pddf.loc[i, "Breakpoint1"] = in_fus_detect_pddf.loc[i, "bp1_lifted"]
                in_fus_detect_pddf.loc[i, "Breakpoint2"] = in_fus_detect_pddf.loc[i, "bp2_lifted"]
        # finally, drop the added columns again and write to disk
        in_fus_detect_pddf.drop(labels=["bp1_lifted", "bp2_lifted"], axis=1, inplace=True)
        in_fus_detect_pddf.to_csv(self.in_fus_detect, sep=";", index=False)
        # urla - note: one could also drop the old cols, rename the new ones and reorder the df, but I find this much more fail-save

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input', dest='input', help='Detected fusions file', required=True)
    parser.add_argument('-c', '--config', dest='config', help='Main config file', required=True)
    parser.add_argument('-l', '--logger', dest='logger', help='Logging of processing steps', default="")
    args = parser.parse_args()

    flo = FusionLiftover(args.input, Config(args.config), args.logger)
    flo.liftcoords()

if __name__ == '__main__':
    main()
