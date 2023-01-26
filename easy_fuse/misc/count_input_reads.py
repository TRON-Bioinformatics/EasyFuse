#!/usr/bin/env python3


from argparse import ArgumentParser
import os
import subprocess


def _get_input_read_count_from_star(infile):
    """Parses a star output log file to get input read counts from the fastq origin"""
    if not os.path.exists(infile):
        return -1
    with open(infile, "r") as inf:
        for line in inf:
            if line.split("|")[0].strip() == "Number of input reads":
                return int(line.split("|")[1].strip())
    return -1


def _get_input_read_count_from_fastq(infile):
    """Parses input FASTQ to get read count"""
    if not os.path.exists(infile):
        return -1
    ps = subprocess.Popen(("zcat", infile), stdout=subprocess.PIPE)
    result = subprocess.check_output(("wc", "-l"), stdin=ps.stdout)
    return int(result) / 4


def get_input_read_count(infile, file_format):
    """Parses input file to get read count"""
    
    if file_format == "star":
        return _get_input_read_count_from_star(infile)
    elif file_format == "fastq":
        return _get_input_read_count_from_fastq(infile)


def run(infile, informat, outfile):
    with open(outfile, "w") as outf:
        outf.write(str(get_input_read_count(infile, informat)))


def add_count_reads_args(parser):
    parser.add_argument("-i", "--input_file", dest="input_file", help="Specify input file")
    parser.add_argument("-f", "--input_format", dest="input_format", help="Specify input file format")
    parser.add_argument("-o", "--output_file", dest="output_file", help="Specify output file")
    parser.set_defaults(func=count_reads_command)


def count_reads_command(args):
    run(args.input_file, args.input_format, args.output_file)
