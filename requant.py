#!/usr/bin/env python
"""Requantification module
"""

__version__ = '0.1'

import sys

from argparse import ArgumentParser

import pysam

class Classification(object):
    
    def __init__(self, bam_file, bed_file, output_file, junc_cutoff):
        self.bam_file = bam_file
        self.bed_file = bed_file
        self.output_file = output_file
        self.junc_cutoff = junc_cutoff
        self.counts = {}
        self.bps = self.get_breakpoints()

    def get_breakpoints(self):
        bps = {}
        with open(self.bed_file) as bed:
            for line in bed:
                elements = line.rstrip().split("\t")
                fgid = elements[0]
                bp_pos = int(elements[2])
                bps.setdefault(fgid, []).append(bp_pos)

        return bps

    def run(self):
        samfile = pysam.AlignmentFile(self.bam_file, 'rb')

        ref_ids = samfile.references[:-1]

        for ref_id in ref_ids:
            reads_on_ref_l = {}
            reads_on_ref_r = {}
            for read in samfile.fetch(ref_id):
                query_name = read.query_name
                start = read.reference_start
                end = read.reference_end
                if read.is_read1:
                    reads_on_ref_l[query_name] = (start, end)
                else:
                    reads_on_ref_r[query_name] = (start, end)

            bp_pos = self.bps[ref_id][0]

            reads_gene_a = 0
            reads_gene_b = 0
            junc_reads = 0
            span_pairs = 0

            for query_name in reads_on_ref_l:
                (start_r1, end_r1) = reads_on_ref_l[query_name]
                (start_r2, end_r2) = reads_on_ref_r[query_name]


                if start_r1 < bp_pos:
                    reads_gene_a += 1

                if start_r2 < bp_pos:
                    reads_gene_a += 1

                if end_r1 > bp_pos:
                    reads_gene_b += 1

                if end_r2 > bp_pos:
                    reads_gene_b += 1

                if (start_r1 <= bp_pos <= end_r2 or
                    start_r2 <= bp_pos <= end_r1):
                    # Increment spanning pair count
                    span_pairs += 1

                if (start_r1 + self.junc_cutoff) <= bp_pos <= (end_r1 - self.junc_cutoff):
                    # increment junction read count
                    junc_reads += 1
                    
                if (start_r2 + self.junc_cutoff) <= bp_pos <= (end_r2 - self.junc_cutoff):                    
                    # increment junction read count
                    junc_reads += 1
                
            self.counts[ref_id] = (reads_gene_a, reads_gene_b, junc_reads, span_pairs)

        print("FGID;MD5;Type;Reads_Gene_A;Reads_Gene_B;Junc_Reads;Span_Pairs")
        for ref_id in sorted(self.counts):
            ref_id_split = ref_id.rsplit("_", 2)
            fgid = ref_id_split[0]
            context_id = ref_id_split[1]
            type = ref_id_split[2]
            (reads_gene_a, reads_gene_b, junc_reads, span_pairs) = self.counts[ref_id]
            print("{};{};{};{};{};{};{}".format(fgid, context_id, type, reads_gene_a, reads_gene_b, junc_reads, span_pairs))


def main():
    parser = ArgumentParser(description='Extracts information on fusion genes')
    parser.add_argument('-i', '--bam', dest='bam', help='Specify input alignment file (BAM).', required=True)
    parser.add_argument('-b', '--bed', dest='bed', help='Specify reference fusion BED file.', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Specify output file with mapping information', default=sys.stdout)
    parser.add_argument('-c', '--junc-cutoff', dest='junc_cutoff', help='Specify cutoff for junction reads. A Read has to overlap the breakpoint with at least <INT> bases.', default=10)

    args = parser.parse_args()

    c = Classification(args.bam, args.bed, args.output, args.junc_cutoff)

    c.run()


if __name__ == '__main__':
    main()
