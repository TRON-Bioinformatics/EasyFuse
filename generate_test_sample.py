#!/usr/bin/env python

from __future__ import print_function

import os
import gzip

from argparse import ArgumentParser

import pysam

def parse_input_table(infile):
    fusion_list = []
    with open(infile) as inf:
        next(inf)
        for line in inf:
            elements = line.rstrip().split(";")
            fusion_list.append(elements[0].split("_"))
    return fusion_list

def parse_ref_genes(infile):
    genes_map = {}
    with open(infile) as inf:
        for line in inf:
            elements = line.rstrip().split("\t")
            chrom = elements[0]
            start = int(elements[1])
            end = int(elements[2])
            tid = elements[3]
            genes_map[tid] = [chrom, start, end]
    
    return genes_map

def parse_gene_symbols(infile):
    symbols_map = {}
    with open(infile) as inf:
        for line in inf:
            elements = line.rstrip().split("\t")
            tid = elements[0]
            gene_symbol = elements[1]
            symbols_map.setdefault(gene_symbol, []).append(tid)

    return symbols_map

class SampleGenerator:
    def __init__(self, input_table, input_bam, ref_bed, gene_symbols_file, working_dir, prefix):
        self.input_table = input_table
        self.input_bam = input_bam
        self.ref_bed = ref_bed
        self.gene_symbols_file = gene_symbols_file
        self.fusion_loci = []
        self.sample_size = 0
        
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)

        self.fq1 = gzip.open(os.path.join(working_dir, prefix + "_R1.fastq.gz"), "wb")
        self.fq2 = gzip.open(os.path.join(working_dir, prefix + "_R2.fastq.gz"), "wb")


    def check_pair(self, r1, r2):
        qname_r1 = r1.query_name
        rname_r1 = r1.reference_name
        start_r1 = r1.reference_start
        end_r1 = r1.reference_end

        qname_r2 = r2.query_name
        rname_r2 = r2.reference_name
        start_r2 = r2.reference_start
        end_r2 = r2.reference_end


        
        valid_r1 = False
        valid_r2 = False

        for (chrom, min_start, max_end) in self.fusion_loci:
            if chrom == rname_r1 or chrom == rname_r2:
                valid_r1 = min_start < start_r1 < max_end or min_start < end_r1 < max_end
                valid_r2 = min_start < start_r2 < max_end or min_start < end_r2 < max_end

                if valid_r1 or valid_r2:
                    break


        if valid_r1 or valid_r2:
            self.sample_size += 2

            seq_r1 = r1.get_forward_sequence()
            seq_r2 = r2.get_forward_sequence()
            qual_r1 = "".join([chr(x+33) for x in r1.get_forward_qualities()])
            qual_r2 = "".join([chr(x+33) for x in r2.get_forward_qualities()])

            # write to FASTQ
            self.fq1.write("@{} 1:N:0:NNNNNN\n{}\n+\n{}\n".format(qname_r1, seq_r1, qual_r1))
            self.fq2.write("@{} 2:N:0:NNNNNN\n{}\n+\n{}\n".format(qname_r2, seq_r2, qual_r2))

#                    print(rname_r1, start_r1, end_r1)
#                    print(rname_r2, start_r2, end_r2)
#                    print()
    

    def run(self):
        samfile = pysam.AlignmentFile(self.input_bam, 'rb')
        genes_map = parse_ref_genes(self.ref_bed)
        #    print(genes_map)
        symbols_map = parse_gene_symbols(self.gene_symbols_file)
        #    print(symbols_map)
        fusion_list = parse_input_table(self.input_table)

        for fusion in fusion_list:
            for gene in fusion:
                min_start = 1000000000
                max_end = 0
                chrom = ""
                for trans in symbols_map[gene]:
                    chrom, start, end = genes_map[trans]
                    if start < min_start:
                        min_start = start
                    if end > max_end:
                        max_end = end
                self.fusion_loci.append((chrom, min_start, max_end))


        print(self.fusion_loci)
        c = 0
        for read in samfile.fetch(until_eof=True):
            flag = read.flag
            if flag > 255:
                continue
        
            c += 1

            if c % 2 == 1:
                r1 = read
            else:
                r2 = read
            
                if flag & 0x40:
                    self.check_pair(r2, r1)
                else:
                    self.check_pair(r1, r2)

            if c % 1000000 == 0:
                print("%s Reads processed" % c)
                print("Current sample size: %s" % self.sample_size)
        
        self.fq1.close()
        self.fq2.close()





def main():
    parser = ArgumentParser(description="Generate Test Sample for EasyFuse pipeline")
    parser.add_argument("-i", "--input-table", dest="input_table", help="Path to input table with fusion gene info", required=True)
    parser.add_argument("-b", "--input-bam", dest="input_bam", help="Path to input BAM file", required=True)
    parser.add_argument("-k", "--reference-genes", dest="ref_genes", help="Path to reference genes table", required=True)
    parser.add_argument("-g", "--gene-symbols", dest="gene_symbols", help="Path to gene symbol mapping table", required=True)
    parser.add_argument("-w", "--working-dir", dest="working_dir", help="Specify working directory for output files", required=True)
    parser.add_argument("-p", "--prefix", dest="prefix", help="Specify output prefix for fastq files", default="test")

    args = parser.parse_args()
    sg = SampleGenerator(args.input_table, args.input_bam, args.ref_genes, args.gene_symbols, args.working_dir, args.prefix)
    sg.run()


if __name__ == "__main__":
    main()
