#!/usr/bin/env python3
from argparse import ArgumentParser


class Parser(object):
    def __init__(self):
        pass

    def parse_total_sequences(self, infile):
        total_sequences = 0
        with open(infile) as inf:
            for line in inf:
                if line.startswith("Total Sequences"):
                    elements = line.rstrip().split("\t")
                    total_sequences = elements[1]
                    return total_sequences
        return total_sequences

    def parse_quality(self, infile, outfile):
        seq_map = {}
        read_length = 0
        with open(infile) as inf:
            check = False
            for line in inf:
                if check and line.startswith(">>END_MODULE"):
                    check = False
                if check and not line.startswith(">>END_MODULE"):
                    elements = line.rstrip().split("\t")
                    base = elements[0]
                    if base != "#Base":
                        mean = float(elements[1])
                        median = float(elements[2])
                        lower_quartile = float(elements[3])  # min 20, so suggested trim length is not too strict
                        upper_quartile = float(elements[4])
                        percentile_10 = float(elements[5])
                        percentile_90 = float(elements[6])
                        if read_length < int(base):
                            read_length = int(base)
                        seq_map[int(base)] = [
                            mean,
                            median,
                            lower_quartile,
                            upper_quartile,
                            percentile_10,
                            percentile_90
                        ]
                if line.startswith(">>Per base sequence quality"):
                    check = True

        counter = 0
        stop = False
        badass_bases = []
        for key in sorted(seq_map.keys(), reverse=True):
            if seq_map[key][2] > 20:
                stop = True
            if not stop:
                counter += 1
            else:
                if seq_map[key][2] < 20:
                    badass_bases.append(key)

        remaining = read_length - counter
        actual = 0
        if (read_length - counter) < (0.75 * read_length):
            actual = int(round(0.25 * read_length))
        else:
            actual = counter
        total_sequences = self.parse_total_sequences(infile)
        # TODO: introduce something better to write a CSV file, eg: pandas or some CSV library
        with open(outfile, 'a') as outf:
            outf.write(
                "{},{},{},{},{},{},{}\n".format(infile, read_length, counter, remaining, actual, len(badass_bases),
                                                total_sequences))
        return seq_map


def add_qc_parser_args(parser):
    parser.add_argument('-i', '--input', dest='input', nargs='+',
                        help='Specify input file(s) (fastqc_data.txt in fastqc folder(s))')
    parser.add_argument('-o', '--qc-table', dest='qc_table', help='Specify output QC table')
    parser.set_defaults(func=qc_parser_command)


def qc_parser_command(args):
    open(args.qc_table, "w").write(
        "filename,read_length,suggested_trim_length,remaining_read_length,actual_trim_length,badass_bases,total_sequences\n")
    p = Parser()

    for file in args.input:
        qmap = p.parse_quality(file, args.qc_table)
