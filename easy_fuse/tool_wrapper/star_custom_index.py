#!/usr/bin/env python


from argparse import ArgumentParser
import math
import subprocess


def get_sequence_count_and_len(fasta_file):
    """Returns the number of sequences and the total bases from a fasta file"""
    seq_num = 0
    total_bases = 0
    with open(fasta_file) as fasta:
        for line in fasta:
            if line[0] == ">":
                seq_num += 1
            else:
                total_bases += len(line.rstrip())
    return seq_num, total_bases


def get_pseudo_genome_adjustments_for_star(fasta_file):
    """Return the genome size of an associated fasta file"""
    (seq_num, genome_size) = get_sequence_count_and_len(fasta_file)
    star_genome_chr_bin_n_bits = min(18, int(math.log(genome_size / seq_num, 2)))
    star_genome_sa_index_n_bases = min(14, int(math.log(genome_size, 2) / 2 - 1)) - 2

    return str(star_genome_chr_bin_n_bits), str(star_genome_sa_index_n_bases)


def run(input_fasta, genome_dir, threads, star_binary):
    (chr_bin_n_bits, sa_index_n_bases) = get_pseudo_genome_adjustments_for_star(input_fasta)
    
    cmd_starindex = "{0} --runMode genomeGenerate --runThreadN {1} --limitGenomeGenerateRAM 48000000000 --genomeChrBinNbits {2} --genomeSAindexNbases {3} --genomeDir {4} --genomeFastaFiles {5}".format(star_binary, threads, chr_bin_n_bits, sa_index_n_bases, genome_dir, input_fasta)

    subprocess.run(cmd_starindex.split(" "))


def add_star_custom_index_args(parser):
    parser.add_argument("-i", "--input_fasta", dest="input_fasta", help="Specify input fasta file", required=True)
    parser.add_argument("-o", "--genome_dir", dest="genome_dir", help="Specify path to output genome dir", required=True)
    parser.add_argument("-t", "--threads", dest="threads", help="Specify amount of threads to be used during runtime", default=1)
    parser.add_argument("-b", "--star_binary", dest="star_binary", help="Specify path to STAR binary for execution", default="STAR")
    parser.set_defaults(func=star_custom_index_command)

def star_custom_index_command(args):
    run(args.input_fasta, args.genome_dir, args.threads, args.star_binary)
