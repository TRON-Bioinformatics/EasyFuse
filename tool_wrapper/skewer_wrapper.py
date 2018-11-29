#!/usr/bin/env python

import sys
import os
import subprocess

from argparse import ArgumentParser

sys.path.append(os.path.join(os.path.dirname( __file__ ), '..'))

from misc.config import Config

def main():
    parser = ArgumentParser(description='Wrapper around skewer.')
    parser.add_argument('-q', '--qc-table', dest='qc_table', help='Specify input QC table', required=True)
    parser.add_argument('-i', '--input', dest='input', nargs='+', help='Specify input FASTQ files', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Specify output folder', default='.')
    parser.add_argument('-m', '--min-read-length-perc', dest='min_read_length_perc', type=float, help='Specify minimum read length percentage', default=0.75)
    args = parser.parse_args()

    cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini'))

    min_rl_perc = args.min_read_length_perc
    remaining_len = 1000
    read_len = 0
    with open(args.qc_table) as inf:
        next(inf)
        for line in inf:
            elements = line.rstrip().split(",")
            read_len = int(elements[1])
            remaining = int(elements[3])
            if remaining < remaining_len:
                remaining_len = remaining

    if remaining_len != read_len and remaining_len >= (min_rl_perc * read_len):
        cmd = "{} -l {} -q 28 -m pe -t 6 -k 5 -z -o out_file {}".format(cfg.get('commands', 'skewer_cmd'), remaining_len, " ".join(args.input))

        proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        (stdoutdata, stderrdata) = proc.communicate()
        r = proc.returncode
        if r != 0:
            sys.stderr.write("Error within skewer\n")
            sys.stderr.write(stderrdata)
            sys.exit(1)
#        else:
#            for line in stdoutdata:
#                sys.stdout.write(line)
#            for line in stderrdata:
#                print line

    elif remaining_len < min_rl_perc * read_len:
        sys.stderr.write("Trim length too long. Discarding fastqs ({}).".format(args.input))
        open(os.path.join(args.output, "to_be_discarded"), 'a').close()
        sys.exit(1)
    elif remaining_len == read_len:
        if len(args.input) > 1:
            out_1 = os.path.join(args.output, "out_file-trimmed-pair1.fastq.gz")
            out_2 = os.path.join(args.output, "out_file-trimmed-pair2.fastq.gz")
            sys.stdout.write("Nothing to trim. Creating symlinks\n")
            sys.stdout.write("ln -s {} {}\n".format(args.input[0], out_1))
            sys.stdout.write("ln -s {} {}\n".format(args.input[1], out_2))
            os.symlink(args.input[0], out_1)
            os.symlink(args.input[1], out_2)
        else:
            out = os.path.join(args.output, "out_file-trimmed.fastq.gz")
            sys.stdout.write("Nothing to trim. Creating symlink\n")
            sys.stdout.write("ln -s {} {}\n".format(args.input[0], out))
            os.symlink(args.input[0], out)


if __name__ == '__main__':
    main()
