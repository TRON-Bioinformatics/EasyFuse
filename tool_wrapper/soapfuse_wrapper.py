#!/usr/bin/env python

import sys
import os
import gzip
import subprocess

from argparse import ArgumentParser


def main():
    parser = ArgumentParser(description="Wrapper around SOAPfuse.")
    parser.add_argument("-q","--qc-table", dest="qc_table", help="Specify input QC table")
    parser.add_argument("-i","--input", dest="input", nargs="+", help="Specify input FASTQ files", required=True)
    parser.add_argument("-o","--output", dest="output", help="Specify output folder", default=".")
    parser.add_argument("-b", "--binary", dest="binary", help="Specify Soapfuse binary", required=True)
    parser.add_argument("-c", "--config-file", dest="config_file", help="Specify Soapfuse config file to use for your analysis", required=True)
    args = parser.parse_args()

    soapfuse_bin = args.binary
    soapfuse_cfg = args.config_file

    cols = args.input[0].rsplit("/", 3)
    remaining_len = 1000
    read_len = 0
    if args.qc_table != "None":

        with open(args.qc_table) as inf:
            inf.readline()
            for line in inf:
                elements = line.rstrip().split(",")
                read_len = int(elements[1])
                remaining = int(elements[3])
                if remaining < remaining_len:
                    remaining_len = remaining
    else:
        with gzip.open(args.input[0]) as inf:
            next(inf)
            line = inf.readline()
            remaining_len = len(line.rstrip())

    
    cfg_file = os.path.join(args.output, "samples.csv")
    cfg_out = open(cfg_file, "w")

    cfg_out.write("{}\t{}\tS\t{}\n".format(cols[1], cols[2], remaining_len))
    cfg_out.close()
    for i, file in enumerate(args.input):
        try:
            os.symlink(file, os.path.join(cols[0], cols[1], cols[2], "S_{}.fastq.gz".format((i+1))))
        except:
            print("Symlink already exists!")
    cmd = "{} -c {} -fd {} -l {} -o {}".format(soapfuse_bin, soapfuse_cfg, cols[0], cfg_file, args.output)

    proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    (stdoutdata, stderrdata) = proc.communicate()
    r = proc.returncode
    if r != 0:
        print(stderrdata)
        sys.exit(1)



if __name__ == "__main__":
    main()
