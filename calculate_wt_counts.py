#!/usr/bin/env python

from __future__ import print_function

import os
import gzip
import time

from argparse import ArgumentParser

import pysam



def parse_input_table(infile):
    fusion_dict = {}
    with open(infile) as inf:
        next(inf)
        for line in inf:
            elements = line.rstrip().split(";")
            fgid = elements[0]
            bp1 = elements[2]
            bp2 = elements[3]
            (chrom1, pos1, strand1) = bp1.split(":")
            (chrom2, pos2, strand2) = bp2.split(":")
            chrom1 = chrom1.strip("chr")
            chrom2 = chrom2.strip("chr")
            fusion_dict[fgid] = (chrom1, pos1, chrom2, pos2)
#            fusion_dict[fgid + "_wt2"] = (chrom2, pos2)
#            fusion_list.append((fgid + "_wt1", chrom1, pos1))
#            fusion_list.append((fgid + "_wt2", chrom2, pos2))
    return fusion_dict


class WTRequantification:
    def __init__(self, input_table, input_bam, output):
        self.input_table = input_table
        self.input_bam = input_bam
        self.fusion_loci = []
        self.sample_size = 0
        self.wt_counts = {}
        self.outf = open(output, 'w')


    def check_pair(self, r1, r2):
        qname_r1 = r1.query_name
        rname_r1 = r1.reference_name
        start_r1 = r1.reference_start
        end_r1 = r1.reference_end

        qname_r2 = r2.query_name
        rname_r2 = r2.reference_name
        start_r2 = r2.reference_start
        end_r2 = r2.reference_end


        for (fgid, chrom1, pos1, chrom2, pos2) in self.fusion_loci:
            if chrom1 == rname_r1:
                if (start_r1+10) <= pos1 <= (end_r1-10):
                    self.wt_counts[fgid][0] += 1
            if chrom1 == rname_r2:
                if (start_r2+10) <= pos1 <= (end_r2-10):
                    self.wt_counts[fgid][0] += 1

            if chrom2 == rname_r1:
                if (start_r1+10) <= pos2 <= (end_r1-10):
                    self.wt_counts[fgid][1] += 1
            if chrom2 == rname_r2:
                if (start_r2+10) <= pos2 <= (end_r2-10):
                    self.wt_counts[fgid][1] += 1
    

    def run(self):
        samfile = pysam.AlignmentFile(self.input_bam, 'rb')

        fusion_list = parse_input_table(self.input_table)

        for fgid in fusion_list:
            self.fusion_loci.append((fgid, fusion_list[fgid][0], int(fusion_list[fgid][1]), fusion_list[fgid][2], int(fusion_list[fgid][3])))
            self.wt_counts[fgid] = [0, 0]

#        print(self.fusion_loci)
        c = 0
        for read in samfile.fetch(until_eof=True):
            flag = read.flag
            if flag > 255:
                continue
        
            c += 1

            if c % 2 == 1:
                r1 = read
            else:
                r2 = read
            
                if flag & 0x40:
                    self.check_pair(r2, r1)
                else:
                    self.check_pair(r1, r2)

            if c % 1000000 == 0:
                print("{} {} Reads processed".format(time.strftime("%Y-%m-%d %H:%M:%S"), c))


        self.outf.write("FGID;wt1_junc_reads;wt2_junc_reads\n")
        for fgid in sorted(self.wt_counts):
            self.outf.write("{};{};{}\n".format(fgid, self.wt_counts[fgid][0], self.wt_counts[fgid][1]))
        self.outf.close()


def main():
    parser = ArgumentParser(description="Generate Test Sample for EasyFuse pipeline")
    parser.add_argument("-i", "--input-table", dest="input_table", help="Path to input table with fusion gene info", required=True)
    parser.add_argument("-b", "--input-bam", dest="input_bam", help="Path to input BAM file", required=True)
    parser.add_argument("-o", "--output", dest="output_table", help="Specify path to output table", required=True)

    args = parser.parse_args()
    sg = WTRequantification(args.input_table, args.input_bam, args.output_table)
    sg.run()


if __name__ == "__main__":
    main()
