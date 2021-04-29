#!/usr/bin/env python

import sys
import os
import gzip
import subprocess

from argparse import ArgumentParser

sys.path.append(os.path.join(os.path.dirname( __file__ ), '..'))

from misc.config import Config

def main():
    parser = ArgumentParser(description='Wrapper around SOAPfuse.')
    parser.add_argument('-q','--qc-table', dest='qc_table', help='Specify input QC table')
    parser.add_argument('-i','--input', dest='input', nargs='+', help='Specify input FASTQ files', required=True)
    parser.add_argument('-o','--output', dest='output', help='Specify output folder', default='.')
    parser.add_argument('-g','--genome', dest='genome', help='Specify reference genome', default='hg38')
    parser.add_argument('-c','--config', dest='config', help='Specify config file')
    args = parser.parse_args()

    if args.config:
        cfg = Config(args.config)
    else:
        cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini'))

    soapfuse_cfg = cfg.get('otherFiles', 'ensembl_soapfuse_cfg_{}'.format(args.genome))

    cols = args.input[0].rsplit("/",3)
#    print(cols)
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
#        print(args.input[0])
        with gzip.open(args.input[0]) as inf:
            next(inf)
            line = inf.readline()
            remaining_len = len(line.rstrip())

#    print(remaining_len)
    
    cfg_file = os.path.join(args.output, "samples.csv")
    cfg_out = open(cfg_file, "w")

    cfg_out.write("{}\t{}\tS\t{}\n".format(cols[1], cols[2], remaining_len))
    cfg_out.close()
    for i, file in enumerate(args.input):
        try:
            os.symlink(file, os.path.join(cols[0], cols[1], cols[2], "S_{}.fastq.gz".format((i+1))))
        except:
            print("Symlink already exists!")
    cmd = "{} -c {} -fd {} -l {} -o {}".format(cfg.get('commands', 'soapfuse_cmd'), soapfuse_cfg, cols[0], cfg_file, args.output)
#    print cmd
    proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    (stdoutdata, stderrdata) = proc.communicate()
    r = proc.returncode
    if r != 0:
        print(stderrdata)
        sys.exit(1)
    else:
        for line in stdoutdata.split("\n"):
            print(line)
        for line in stderrdata.split("\n"):
            print(line)



if __name__ == '__main__':
    main()
