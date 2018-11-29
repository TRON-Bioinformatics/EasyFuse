#!/usr/bin/env python

import sys
import os
import subprocess

from argparse import ArgumentParser

sys.path.append(os.path.join(os.path.dirname( __file__ ), '..'))

from misc.config import Config

def main():
    parser = ArgumentParser(description='Wrapper around MapSplice.')
    parser.add_argument('-i', '--input', dest='input', nargs='+', help='Specify input FASTQ files', required=True)
    parser.add_argument('-t', '--threads', dest='threads', help='Specify number of threads to be used')
    parser.add_argument('-o', '--output', dest='output', help='Specify output folder', default='.')
    parser.add_argument('-g', '--reference-genome', dest='genome', help='Specify reference genome', default='hg38')
    parser.add_argument('-c', '--config', dest='config', help='Specify config file')
    args = parser.parse_args()

    if args.config:
        cfg = Config(args.config)
    else:
        cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini'))
    
    
    tmp_file_1 = os.path.join(args.output, "S_64phred_R1_001.fastq")
    tmp_file_2 = os.path.join(args.output, "S_64phred_R2_001.fastq")
    cmd_1 = "{} fastq-format --phred33to64 --strip --suffix /1 -in <(zcat {}) --out {} > {}".format(cfg.get('commands','phred33to64_cmd'), args.input[0], tmp_file_1, os.path.join(args.output,"mapsplice_prep1.log"))
    cmd_2 = "{} fastq-format --phred33to64 --strip --suffix /2 -in <(zcat {}) --out {} > {}".format(cfg.get('commands','phred33to64_cmd'), args.input[1], tmp_file_2, os.path.join(args.output,"mapsplice_prep2.log"))

    proc = subprocess.Popen(['/bin/bash', '-c', cmd_1], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    (stdoutdata, stderrdata) = proc.communicate()
    r = proc.returncode
    if r != 0:
        print("Error within phred64 conversion (R1)")
        print(stderrdata)
        sys.exit(1)

    proc = subprocess.Popen(['/bin/bash', '-c', cmd_2], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    (stdoutdata, stderrdata) = proc.communicate()
    r = proc.returncode
    if r != 0:
        print("Error within phred64 conversion (R2)")
        print(stderrdata)
        sys.exit(1)

    cmd = "{} -c {} -x {} -1 {} -2 {} -p {} -o {} --gene-gtf {} --fusion".format(cfg.get('commands', 'mapsplice_cmd'), cfg.get('bowtie_indexes', 'bowtie_index_{}_fasta'.format(args.genome)), cfg.get('bowtie_indexes', 'bowtie_index_{}'.format(args.genome)), tmp_file_1, tmp_file_2, args.threads, args.output, cfg.get('references', 'ensembl_genes_gtf_{}'.format(args.genome)))

    proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    (stdoutdata, stderrdata) = proc.communicate()
    
    r = proc.returncode
    if r != 0:
        print("Error within mapsplice")
        print(stderrdata)
        sys.exit(1)
    else:
        for line in stdoutdata.split("\n"):
            print(line)
        for line in stderrdata.split("\n"):
            print(line)


if __name__ == '__main__':
    main()
