#!/usr/bin/env python

import sys
from argparse import ArgumentParser

class Parser(object):
    def __init__(self):
        pass

    def parse_total_sequences(self,infile):
        total_sequences = 0
        counter = 0
        with open(infile) as inf:
            for line in inf:
                if line.startswith("Total Sequences"):
                    elements = line.rstrip().split("\t")
                    total_sequences = elements[1]
                    return total_sequences
        return total_sequences

    def parse_overrepresented(self,infile,outfile):
        seq_map = {}
        with open(infile) as inf:
            check = False
            for line in inf:
                if check and line.startswith(">>END_MODULE"):
                    check = False
                if check and not line.startswith(">>END_MODULE"):
#                    print line
                    elements = line.rstrip().split("\t")
                    sequence = elements[0]
                    occurences = elements[1]
                    perc = elements[2]
                    if sequence != "#Sequence":
                        seq_map[sequence] = ""
                if line.startswith(">>Overrepresented sequences"):
                    check = True
        outf = open(outfile,'a')
        for key in seq_map:
            outf.write(key + "\n")
        outf.close()
        return seq_map

    def parse_quality(self,infile,outfile):
        seq_map = {}
        read_length = 0
        with open(infile) as inf:
            check = False
            for line in inf:
                if check and line.startswith(">>END_MODULE"):
                    check = False
                if check and not line.startswith(">>END_MODULE"):
#                    print line
                    elements = line.rstrip().split("\t")
                    base = elements[0]
                    if base != "#Base":
                        mean = float(elements[1])
                        median = float(elements[2])
                        lower_quartile = float(elements[3])       # min 20, so suggested trim length is not too strict
                        upper_quartile = float(elements[4])
                        percentile_10 = float(elements[5])
                        percentile_90 = float(elements[6])
                        if read_length < int(base):
                            read_length = int(base)
                        seq_map[int(base)] = [mean,median,lower_quartile,upper_quartile,percentile_10,percentile_90]
                if line.startswith(">>Per base sequence quality"):
                    check = True

        counter = 0
        stop = False
        badass_bases = []
        for key in sorted(seq_map.keys(),reverse=True):
            if seq_map[key][2] > 20:
                stop = True
            if not stop:
                counter += 1
            else:
                if seq_map[key][2] < 20:
                    badass_bases.append(key)
                    
        remaining = read_length-counter
        actual = 0
#        print "Suggested Trim cutoff: " + str(counter) + " bp"
#        print "Remaining Read length: " + str(remaining) + " bp"
        if (read_length-counter) < (0.75*read_length):
            actual = int(round(0.25*read_length))
#            print "Actual Trim length: " + str(actual)
        else:
            actual = counter
#            print "Actual Trim length: " + str(actual)
#        print "Badass bases: " + str(badass_bases)
        total_sequences = self.parse_total_sequences(infile)
        outf = open(outfile,'a')

        outf.write(infile + "," + str(read_length) + "," + str(counter) + "," + str(remaining) + "," + str(actual) + "," + str(len(badass_bases)) + "," + str(total_sequences) + "\n")
#        print infile + "," + str(counter) + "," + str(remaining) + "," + str(actual) + "," + str(len(badass_bases))
        outf.close()
        return seq_map
 
    def write_fasta(self,outfile,seq_map):
        outf = open(outfile,"w")
        for i,key in enumerate(seq_map,start=1):
            header = ">seq" + str(i)
            sequence = key
            outf.write(header+"\n")
            outf.write(sequence+"\n")

def main():
    parser = ArgumentParser(description='handles garbage in racoon-style')
    parser.add_argument('-i','--input',dest='input',nargs='+',help='Specify input file(s) (fastqc_data.txt in fastqc folder(s))')
    parser.add_argument('-o','--overrepresented',dest='overrepresented',help='Specify output overrepresented sequences file')
    parser.add_argument('-q','--qc-table',dest='qc_table',help='Specify output QC table')
    args = parser.parse_args()
#    infile = sys.argv[1]
#    outfile = sys.argv[2]
    open(args.qc_table,"w").write("filename,read_length,suggested_trim_length,remaining_read_length,actual_trim_length,badass_bases,total_sequences\n")
    p = Parser()

    for file in args.input:
        seq_map = p.parse_overrepresented(file,args.overrepresented)
#        p.write_fasta(outfile,seq_map)
#        total_sequences = p.parse_total_sequences(file)
#        print total_sequences
        qmap = p.parse_quality(file,args.qc_table)
#    print qmap


if __name__ == '__main__':
    main()
