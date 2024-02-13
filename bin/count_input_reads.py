import os
import subprocess


def __get_input_read_count_from_star(infile):
    """Parses a star output log file to get input read counts from the fastq origin"""
    if not os.path.exists(infile):
        return -1
    with open(infile, "r") as inf:
        for line in inf:
            if line.split("|")[0].strip() == "Number of input reads":
                return int(line.split("|")[1].strip())
    return -1


def __get_input_read_count_from_fastq(infile):
    """Parses input FASTQ to get read count"""
    if not os.path.exists(infile):
        return -1
    ps = subprocess.Popen(("zcat", infile), stdout=subprocess.PIPE)
    result = subprocess.check_output(("wc", "-l"), stdin=ps.stdout)
    return int(result) / 4


def get_input_read_count(infile, file_format):
    """Parses input file to get read count"""

    if file_format == "star":
        return __get_input_read_count_from_star(infile)
    elif file_format == "fastq":
        return __get_input_read_count_from_fastq(infile)
