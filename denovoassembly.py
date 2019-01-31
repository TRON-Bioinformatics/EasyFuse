#!/usr/bin/env python

"""
Identify and extract reads belonging to fusion events
Run a de novo assembly on those reads with trinity
Blast assembled contigs against human transcripts

@author: BNT (URLA)
@version: 20181126
"""

from __future__ import print_function
from argparse import ArgumentParser
import os.path
import sys
import time
import shutil
import shlex
import pysam # pysam is not available for windows (where I run pylint) => pylint: disable=E0401
from easyfuse_scheduling import Queue
from easyfuse_config import Config
from easyfuse_iomethods import IOmethods as urlaio

# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class Assembler(object):
    """Select alignments belonging to putative fusions from an s/bam file"""
    def __init__(self, input_files, cfg, in_bam, output_dir):
        self.input_files = input_files
        self.cfg = cfg
        self.in_bam = in_bam
        self.output_dir = output_dir
        self.glist = set()
        self.glist_pairs = set()

    def get_gene_list(self, glist_file):
        """Read a file of fusion gene ids and return a list of unique genes"""
        with open(glist_file, "r") as in_list_read:
            for line in in_list_read:
                if "--" in line:
                    self.glist.update(line.rstrip("\n").split("--"))
                    self.glist_pairs.add(line.rstrip("\n"))
                else:
                    self.glist.add(line.rstrip("\n"))

    def get_reads_from_bam(self, gcoords_dict, out_dir):
        """For a list of gene ids + their genomic coordinates, fetch all reads mapping to this gene and write them to a bam"""
        count_lines = 0
        for gene in gcoords_dict:
            count_lines += 1
            samfile = pysam.AlignmentFile(self.in_bam, "rb")
            output = os.path.join(out_dir, "Alignments_to_{}.bam".format(gene))
            filtered_bam = pysam.AlignmentFile(output, "wb", template=samfile)
            for read in samfile.fetch(gcoords_dict[gene][0], gcoords_dict[gene][1], gcoords_dict[gene][2]):
                filtered_bam.write(read)

    def get_coords_from_parsedgtf(self, gcoords_file, gcoords_dict):
        """Get genomic coordinates of a gene from a pre-parsed gtf file (parsed with reformatgtf.py)"""
        with open(gcoords_file, "r") as in_gtf_read:
            for line in in_gtf_read:
                line_splitter = line.rstrip("\n").split("\t")
                if line_splitter[0] in self.glist:
                    gcoords_dict[line_splitter[0]] = line_splitter[1], int(line_splitter[2]), int(line_splitter[3])
        return gcoords_dict

    def parse_blast_output(self, in_blast, out_blast, in_fasta, in_fusion):
        """For a gene fusion pair, parse the output of their format-6 custom blast out table and return the fusion transcript and the breakpoint position within it"""
        fusion_pair = in_fusion.split("--")
        blast_dict = {}
        count_fp_hits = 0
        with open(in_blast, "r") as blastf:
            for line in blastf:
                line_splitter = line.rstrip("\n").split("\t")
                #query: TRINITY_DN17_c4_g3_i1; len BRAF: 4; len AKAP9: 7
                try:
                    fus_idx = fusion_pair.index(line_splitter[1].split("_")[0])
                    if not line_splitter[0] in blast_dict:
                        line_buffer = []
                        blast_dict[line_splitter[0]] = [0, 0, line_buffer]
                    blast_dict[line_splitter[0]][fus_idx] += 1
                    blast_dict[line_splitter[0]][2].append(line.rstrip("\n").split("\t") + [fus_idx])
                except ValueError:
                    count_fp_hits += 1

        #print(blast_dict)
        for query in blast_dict:
            if blast_dict[query][0] > 0 and blast_dict[query][1] > 0:
                fus_head, fus_seq = self.get_fasta_record(query, in_fasta)
#                fus1_dict = {}
#                fus2_dict = {}
                breakpoint_pos = 0
                for blast_line_list1 in blast_dict[query][2]:
                    if blast_line_list1[-1] == 1:
                        continue
                    fus1_start = int(blast_line_list1[9])
                    fus1_stop = int(blast_line_list1[10])
                    for blast_line_list2 in blast_dict[query][2]:
                        if blast_line_list2[-1] == 0:
                            continue
                        fus2_start = int(blast_line_list2[9])
                        fus2_stop = int(blast_line_list2[10])
                        if fus1_stop + 1 == fus2_start:
                            if breakpoint_pos > 0 and breakpoint_pos != fus1_stop:
                                print("Warning: Multiple breakpoints detected! Though this might be real, it is more likely to be a technical artifact!")
                            else:
                                breakpoint_pos = fus1_stop
                            #print("--".join([query, str(fus1_stop), blast_line_list1[1], blast_line_list2[1]]))
                        elif fus2_stop + 1 == fus1_start:
                            if breakpoint_pos > 0 and breakpoint_pos != fus2_stop:
                                print("Warning: Multiple breakpoints detected! Though this might be real, it is more likely to be a technical artifact!")
                            else:
                                breakpoint_pos = fus2_stop
                            #print("--".join([query, str(fus2_stop), blast_line_list2[1], blast_line_list1[1]]))
                with open(out_blast, "a") as outf:
                    outf.write(fus_head + "\n")
                    outf.write(fus_seq + "\n")
                    outf.write(fus_seq[breakpoint_pos-100:breakpoint_pos] + "\n")
                    outf.write(fus_seq[breakpoint_pos:breakpoint_pos+100] + "\n")

    @staticmethod
    def get_fasta_record(seq_id, in_fasta):
        """Return the fasta record of a defined sequence id as a tuple (header, sequence)"""
        #in_fasta = "/mnt/bfx/urla/PRJNA252360_syntheticSpikeInFusions/fusionProcessing100bp_filtered/Sample_SRR1659951/fetchdata/fd_1_tool/assembly/trinity_AKAP9_BRAF/Trinity.fasta"
        out_fasta = ["", ""]
        found_rec = False
        with open(in_fasta, "r") as fasf:
            for fasta_line in fasf:
                if fasta_line.startswith(">"):
                    if seq_id in fasta_line:
                        out_fasta[0] = fasta_line[1:].rstrip("\n")
                        found_rec = True
                    else:
                        found_rec = False
                elif found_rec:
                    out_fasta[1] += fasta_line.rstrip("\n")
        return out_fasta

    @staticmethod
    def concatenate_gzip_files(gzip1, gzip2, gzip12):
        """Concatenate two gzip files"""
        with open(gzip12, 'wb') as gf12:
            with open(gzip1, 'rb') as gf1:
                shutil.copyfileobj(gf1, gf12)
            with open(gzip2, 'rb') as gf2:
                shutil.copyfileobj(gf2, gf12)

    def run(self, ref_trans, ref_genome):
        """Process the whole assembly pipeline"""

        # Creating a set of ids of genes involved in a fusion
        start_time = time.time()
        print("Reading gene ids from provided input files...")
        for in_files in self.input_files:
            self.get_gene_list(in_files)
        print("done. Finished in {} seconds.".format(time.time() - start_time))

        # Creating a dictionary of gene ids as keys and their genomic coordinates as list like [chr, start, stop]
        start_time = time.time()
        print("Retrieving genomic coordinates for the input genes...")
        gcoords_file = self.cfg.get("references", "{0}_genes_gtf_{1}".format(ref_trans, ref_genome)) + ".togenecoords"
        if not os.path.isfile(gcoords_file):
            print("file not found")
            sys.exit(1)
        gcoords_dict = {}
        gcoords_dict = self.get_coords_from_parsedgtf(gcoords_file, gcoords_dict)
        print("done. Finished in {} seconds.".format(time.time() - start_time))

        # Fetch reads mapping to a fusion gene from a bam with pysam
        start_time = time.time()
        print("Fetching reads from the original bam file...")
        bam_folder = os.path.join(self.output_dir, "bam_files")
        read_folder = os.path.join(self.output_dir, "read_files")
        urlaio.create_folder_wo_logging(bam_folder)
        urlaio.create_folder_wo_logging(read_folder)

        # check if a coordinate sorted version of the bam file already exists
        if os.path.isfile("{}.sorted.bam".format(self.in_bam[:-4])):
            print("Looks like a coordinate sorted file is already present. Going to use this as input...")
            self.in_bam = "{}.sorted.bam".format(self.in_bam[:-4])

        if "sorted" not in self.in_bam:
            print("The input bam doesn't seem to be sorted. Sorting and indexing now...")
            cmd_sort = "{0} sort -o {1}.sorted.bam {1}.bam".format(self.cfg.get('commands', 'samtools_cmd'), self.in_bam[:-4])
            cmd_index = "{0} index {1}.sorted.bam".format(self.cfg.get('commands', 'samtools_cmd'), self.in_bam[:-4])
            Queue.submit_nonqueue(cmd_sort.split(" "))
            Queue.submit_nonqueue(cmd_index.split(" "))
            self.in_bam = "{}.sorted.bam".format(self.in_bam[:-4])
            print("done. Finished sorting and indexing in {} seconds.".format(time.time() - start_time))
            start_time = time.time()

        self.get_reads_from_bam(gcoords_dict, bam_folder)
        print("done. Finished in {} seconds.".format(time.time() - start_time))

        # Running samtools to extract reads from bam files
        start_time = time.time()
        print("Processing individual bam files and preparing trinity input...")
        print("Extracting reads from bam files...")
        for gene in self.glist:
            gene_bam = os.path.join(bam_folder, "Alignments_to_{}.bam".format(gene))
            gene_read0 = os.path.join(read_folder, "Reads_to_{}_R0.fastq.gz".format(gene))
            gene_read1 = os.path.join(read_folder, "Reads_to_{}_R1.fastq.gz".format(gene))
            gene_read2 = os.path.join(read_folder, "Reads_to_{}_R2.fastq.gz".format(gene))
            cmd_nsort = "{0} sort -n -o {1}.nSorted.bam {1}.bam".format(self.cfg.get('commands', 'samtools_cmd'), gene_bam[:-4])
            cmd_fixmate = "{0} fixmate {1}.nSorted.bam {1}_fixmate.bam".format(self.cfg.get('commands', 'samtools_cmd'), gene_bam[:-4])
            cmd_bam_to_fastq = "{0} fastq -N -0 {1} -1 {2} -2 {3} --threads 4 {4}_fixmate.bam".format(self.cfg.get('commands', 'samtools_cmd'), gene_read0, gene_read1, gene_read2, gene_bam[:-4])

            Queue.submit_nonqueue(cmd_nsort.split(" "))
            Queue.submit_nonqueue(cmd_fixmate.split(" "))
            Queue.submit_nonqueue(cmd_bam_to_fastq.split(" "))

        print("Concatenating reads mapping to fusion partners...")
        # Running "cat" to concatenate read files from each of the fusion partner
        for gene_pair in self.glist_pairs:
            gene1, gene2 = gene_pair.split("--")
            gene1_read1 = os.path.join(read_folder, "Reads_to_{}_R1.fastq.gz".format(gene1))
            gene2_read1 = os.path.join(read_folder, "Reads_to_{}_R1.fastq.gz".format(gene2))
            gene12_read1 = os.path.join(read_folder, "Reads_to_{}_{}_R1.fastq.gz".format(gene1, gene2))
            self.concatenate_gzip_files(gene1_read1, gene2_read1, gene12_read1)
            #cmd_concat_fastq1 = "cat {0} {1} > {2}".format(gene1_read1, gene2_read1, gene12_read1)

            gene1_read2 = os.path.join(read_folder, "Reads_to_{}_R2.fastq.gz".format(gene1))
            gene2_read2 = os.path.join(read_folder, "Reads_to_{}_R2.fastq.gz".format(gene2))
            gene12_read2 = os.path.join(read_folder, "Reads_to_{}_{}_R2.fastq.gz".format(gene1, gene2))
            self.concatenate_gzip_files(gene1_read2, gene2_read2, gene12_read2)
            #cmd_concat_fastq2 = "cat {0} {1} > {2}".format(gene1_read2, gene2_read2, gene12_read2)

            #Queue.submit_nonqueue(cmd_concat_fastq1.split(" "))
            #Queue.submit_nonqueue(cmd_concat_fastq2.split(" "))
        print("done. Finished in {} seconds.".format(time.time() - start_time))

        # Running trinity for de novo assembly of reads mapping to the fusion partners
        start_time = time.time()
        print("Running trinity on fusion reads...")
        trinity_root_folder = os.path.join(self.output_dir, "trinity")
        urlaio.create_folder_wo_logging(trinity_root_folder)
        for gene_pair in self.glist_pairs:
            trinity_gene_pair_folder = os.path.join(trinity_root_folder, "trinity_" + gene_pair.replace("--", "_"))
            urlaio.create_folder_wo_logging(trinity_gene_pair_folder)
            gene1, gene2 = gene_pair.split("--")
            gene12_read1 = os.path.join(read_folder, "Reads_to_{}_{}_R1.fastq.gz".format(gene1, gene2))
            gene12_read2 = os.path.join(read_folder, "Reads_to_{}_{}_R2.fastq.gz".format(gene1, gene2))
            cmd_trinity = "{0} --seqType fq --max_memory 20G --left {1} --right {2} --CPU 8 --output {3}".format(self.cfg.get('commands', 'trinity_cmd'), gene12_read1, gene12_read2, trinity_gene_pair_folder)

            Queue.submit_nonqueue(cmd_trinity.split(" "))
        print("done. Finished in {} seconds.".format(time.time() - start_time))

        # Blasting de novo assembled transcripts against all human transcripts
        start_time = time.time()
        print("Blasting transcripts from trinity against all human transcripts...")
        blast_db = self.cfg.get("references", "{0}_genes_fasta_{1}".format(ref_trans, ref_genome)) + ".toNameHeader"
        blast_folder = os.path.join(self.output_dir, "blast")
        urlaio.create_folder_wo_logging(blast_folder)
        for gene_pair in self.glist_pairs:
            trinity_transcripts = os.path.join(trinity_root_folder, "trinity_" + gene_pair.replace("--", "_"), "Trinity.fasta")
            blast_output = os.path.join(blast_folder, gene_pair.replace("--", "_") + ".blast.out")
            cmd_blast = "{0} -db {1} -query {2} -outfmt \"6 qseqid sseqid pident qlen slen length mismatch gapope evalue bitscore qstart qend sstart send\" -out {3}".format(self.cfg.get('commands', 'blastn_cmd'), blast_db, trinity_transcripts, blast_output)

            Queue.submit_nonqueue(shlex.split(cmd_blast))
        print("done. Finished in {} seconds.".format(time.time() - start_time))

        # Processing blast results and parsing into final data output
        start_time = time.time()
        print("Parsing blast output and extracting final fusion transcripts...")
        for gene_pair in self.glist_pairs:
            in_blast = os.path.join(blast_folder, gene_pair.replace("--", "_") + ".blast.out")
            out_blast = os.path.join(blast_folder, gene_pair.replace("--", "_") + ".blast.out.fusionTranscript")
            in_fasta = os.path.join(trinity_root_folder, "trinity_" + gene_pair.replace("--", "_"), "Trinity.fasta")
            self.parse_blast_output(in_blast, out_blast, in_fasta, gene_pair)
        print("done. Finished in {} seconds.".format(time.time() - start_time))


def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input', dest='input', nargs='+', help='Input files containing fusion genes', required=True)
    parser.add_argument('-c', '--config', dest='config', help='Specify config file.', required=True)
    parser.add_argument('-b', '--bam', dest='bam', help='BAM file of initial read mappings.', required=True)
    parser.add_argument('-t', '--reference-transcripts', dest='ref_trans', choices=['ucsc', 'refseq', 'ensembl'], help='Specify the reference transcript database name.', required=True)
    parser.add_argument('-g', '--reference-genome', dest='ref_genome', choices=['hg18', 'hg19', 'hg38', 'mm9', 'mm10'], help='Specify the reference genome name.', required=True)
    parser.add_argument('-o', '--outdir', dest='outdir', help='Output directory')
    args = parser.parse_args()

    cfg = Config(args.config)
    ass = Assembler(args.input, cfg, args.bam, args.outdir)
    ass.run(args.ref_trans, args.ref_genome)

if __name__ == '__main__':
    main()
