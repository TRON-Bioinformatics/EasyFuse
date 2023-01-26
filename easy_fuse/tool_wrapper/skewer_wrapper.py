#!/usr/bin/env python3

import os
import subprocess
import sys
from shutil import copyfile
from logzero import logger


def run_skewer(skewer_bin, min_rl_perc, qc_table, input, output):

    remaining_len = 1000
    read_len = 0
    with open(qc_table) as inf:
        next(inf)
        for line in inf:
            elements = line.rstrip().split(",")
            read_len = int(elements[1])
            remaining = int(elements[3])
            if remaining < remaining_len:
                remaining_len = remaining

    if remaining_len != read_len and remaining_len >= (min_rl_perc * read_len):
        cmd = "{} -l {} -q 28 -m pe -t 6 -k 5 -z -o {} {}".format(skewer_bin, remaining_len,
                                                                  os.path.join(args.output, "out_file"),
                                                                  " ".join(args.input))

        proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        (stdoutdata, stderrdata) = proc.communicate()
        r = proc.returncode
        if r != 0:
            sys.stderr.write("Error within skewer\n")
            sys.stderr.write(stderrdata)
            sys.exit(1)

    elif remaining_len < min_rl_perc * read_len:
        sys.stderr.write("Trim length too long. Discarding fastqs ({}).".format(args.input))
        open(os.path.join(output, "to_be_discarded"), "a").close()
        sys.exit(1)
    elif remaining_len == read_len:
        if len(input) == 2:
            out_1 = os.path.join(output, "out_file-trimmed-pair1.fastq.gz")
            out_2 = os.path.join(output, "out_file-trimmed-pair2.fastq.gz")
            sys.stdout.write("Nothing to trim. Copying FASTQs\n")
            sys.stdout.write("cp {} {}\n".format(input[0], out_1))
            sys.stdout.write("cp {} {}\n".format(input[1], out_2))
            try:
                copyfile(input[0], out_1)
            #                os.symlink(args.input[0], out_1)
            except OSError:
                print("Copy of {} already exists".format(out_1))
            try:
                copyfile(input[1], out_2)
            #                os.symlink(args.input[1], out_2)
            except OSError:
                logger.error("Copy of {} already exists".format(out_2))


def add_skewer_args(parser):
    parser.add_argument("-q", "--qc-table", dest="qc_table", help="Specify input QC table", required=True)
    parser.add_argument("-i", "--input", dest="input", nargs="+", help="Specify input FASTQ files", required=True)
    parser.add_argument("-o", "--output", dest="output", help="Specify output folder", default=".")
    parser.add_argument("-b", "--binary", dest="binary", help="Specify Skewer binary", required=True)
    parser.add_argument("-m", "--min-read-len-perc", dest="min_read_len_perc", type=float,
                        help="Specify minimum read length percentage", default=0.75)
    parser.set_defaults(func=skewer_command)


def skewer_command(args):
    run_skewer(
        skewer_bin=args.binary,
        min_rl_perc=args.min_read_len_perc,
        qc_table=args.qc_table,
        input=args.input,
        output=args.output)
