#!/usr/bin/env python

import gzip
import os
import subprocess
import sys

from logzero import logger


def run_soapfuse(binary, config_file, qc_table, input, output):

    cols = input[0].rsplit("/", 3)
    remaining_len = 1000
    if qc_table != "None":

        with open(qc_table) as inf:
            inf.readline()
            for line in inf:
                elements = line.rstrip().split(",")
                remaining = int(elements[3])
                if remaining < remaining_len:
                    remaining_len = remaining
    else:
        with gzip.open(input[0]) as inf:
            next(inf)
            line = inf.readline()
            remaining_len = len(line.rstrip())

    cfg_file = os.path.join(output, "samples.csv")
    cfg_out = open(cfg_file, "w")

    cfg_out.write("{}\t{}\tS\t{}\n".format(cols[1], cols[2], remaining_len))
    cfg_out.close()
    for i, file in enumerate(input):
        try:
            os.symlink(file, os.path.join(cols[0], cols[1], cols[2], "S_{}.fastq.gz".format((i + 1))))
        except:
            logger.info("Symlink already exists!")
    cmd = "{} -c {} -fd {} -l {} -o {}".format(binary, config_file, cols[0], cfg_file, output)

    proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    (stdoutdata, stderrdata) = proc.communicate()
    r = proc.returncode
    if r != 0:
        logger.error(stderrdata)
        sys.exit(1)


def add_soapfuse_wrapper_args(parser):
    parser.add_argument("-q", "--qc-table", dest="qc_table", help="Specify input QC table")
    parser.add_argument("-i", "--input", dest="input", nargs="+", help="Specify input FASTQ files", required=True)
    parser.add_argument("-o", "--output", dest="output", help="Specify output folder", default=".")
    parser.add_argument("-b", "--binary", dest="binary", help="Specify Soapfuse binary", required=True)
    parser.add_argument("-c", "--config-file", dest="config_file",
                        help="Specify Soapfuse config file to use for your analysis", required=True)
    parser.set_defaults(func=soapfuse_wrapper_command)


def soapfuse_wrapper_command(args):
    run_soapfuse(
        binary=args.binary,
        config_file=args.config_file,
        qc_table=args.qc_table,
        input=args.input,
        output=args.output
    )
