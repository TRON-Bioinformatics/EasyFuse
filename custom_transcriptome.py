#!/usr/bin/env python

import os
import re
import sys
import time
import csv
import math

from argparse import ArgumentParser

from misc.config import Config
from misc.samples import Samples
import misc.io_methods as IOMethods
import misc.queue as Queueing



def create_folder(folder):
    '''This function creates a specified folder, if not already existing and grants the respective permission.'''
    if not os.path.exists(folder):
        print "Creating folder " + folder
        os.makedirs(folder)
        grant_permissions(folder, stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)

def grant_permissions(path, mode):
    '''This function grants specific permissions for a given folder or file.'''
    for root, dirs, files in os.walk(path, topdown=False):
        for dir in dirs:
            os.chmod(dir, mode)
        for file in files:
            os.chmod(file, mode)

class CustomTranscriptome(object):
    def __init__(self,cfg,input,fusions_table,working_dir,partitions):
#        self.bam_list = []
#        self.job_list = []
        csv.field_size_limit(999999999)
        self.cfg = cfg
        self.input = []
        for i in input:
            self.input.append(i.rstrip("/"))
        self.fusions_table = fusions_table
#        self.min_rl_perc = min_rl_perc
        self.working_dir = working_dir.rstrip("/")
        create_folder(self.working_dir)
        self.samples = Samples(os.path.join(self.working_dir,"samples.csv"))
        self.fasta = os.path.join(self.working_dir,"context.fa")
        self.bed = os.path.join(self.working_dir,"context.bed")
#        self.no_qc = no_qc
        self.partitions = partitions
        self.csv_to_fasta()
        self.generate_bed()
        self.build_index()



    def get_jobs_with_name(self,name):
        '''This function returns a list of job ids for a given sample name.'''
        jobs = Jobs()
        return jobs.get_jobs_with_name_slurm(name)

    def get_fastq_files(self,paths):
        '''This function returns a list of fastq files for the given list of paths.'''
        fastqs = []
        for path in paths:
            if os.path.isdir(path):
                files = os.listdir(path)
                for file in files:
                    file_path = os.path.join(path,file)
                    if os.path.isfile(file_path) and file.endswith((".fastq.gz",".fq.gz",".fq",".fastq")):
                        fastqs.append(file_path)
                    elif os.path.isdir(file_path):
                        fastqs_tmp = self.get_fastq_files([file_path])
                        fastqs.extend(fastqs_tmp)
            elif os.path.isfile(path) and path.endswith((".fastq.gz",".fq.gz",".fq",".fastq")):
                fastqs.append(path)
        return fastqs

    def calc_genome_size(self):
        genome_size = 0
        with open(self.fasta) as inf:
            for line in inf:
                if not line.rstrip().startswith(">"):
                    genome_size += len(line.rstrip())
        print "Custom Genome Size: " + str(genome_size) + "bp"
        return genome_size

### FGID;Fusion_Gene;Breakpoint1;Breakpoint2;FTID;context_sequence_id;type;exon_boundary1;exon_boundary2;exon_boundary;bp1_frame;bp2_frame;frame;wt1_expression;wt2_expression;context_sequence;context_sequence_bp;wt1_context_sequence;wt1_context_sequence_bp;wt2_context_sequence;wt2_context_sequence_bp;neo_peptide_sequence;neo_peptide_sequence_bp;transcript_sequence ###

    def csv_to_fasta(self):
        context_seqs = {}
#        outf = open(self.fasta,"w")
        print self.fusions_table
        with open(self.fusions_table) as csvfile:
            inf = csv.reader(csvfile, delimiter=';')
            header = inf.next()
            print header
#            for row in inf:
#            next(inf)
            col = {}
            for i,colname in enumerate(header):
                col[colname] = i
            for c,row in enumerate(inf):
#                print row
#                elements = line.rstrip().split(";")
                ftid = row[col["FTID"]]
                fgid = row[col["FGID"]]
                context_seq = row[col["context_sequence"]].upper()
                md5_hash = row[col["context_sequence_id"]]
                id = ftid + "_fg_" + str(c)
                context_id = fgid + "_" + md5_hash


#                breakpoint1 = row[col["Breakpoint1"]]
#                breakpoint2 = row[col["Breakpoint2"]]
#                fusion_gene = elements[1]
#                transcript_seq = row[col["transcript_sequence"]].rstrip("NA").upper()
#                rel_bp_1 = row[col["relative_breakpoint1"]]
#                rel_bp_2 = row[col["relative_breakpoint2"]]
#                transcript_sequence = elements[37]
                wt1_context_seq = row[col["wt1_context_sequence"]].rstrip("NA").upper()
                wt2_context_seq = row[col["wt2_context_sequence"]].rstrip("NA").upper()

                context_seqs[context_id] = [context_seq,wt1_context_seq,wt2_context_seq]
#                outf.write(">" + fusion_gene + "_" + breakpoint1 + "_" + breakpoint2 + "_" + str(c) + "\n")
#                for i in range(0,len(context_sequence),60):
#                    outf.write(sequence[i:i+60]+"\n")
#                outf.write(">" + id + "\n")
#                for i in range(0,len(transcript_seq),60):
#                    outf.write(transcript_seq[i:i+60]+"\n")
#        outf.close()

        outf_context = open(self.fasta,"w")
        genome_size = 0
        for context_id in context_seqs:
            (context_seq,wt1_context_seq,wt2_context_seq) = context_seqs[context_id]
            genome_size += (len(context_seq) + len(wt1_context_seq) + len(wt2_context_seq))
            # ft
            outf_context.write(">" + context_id + "_ft\n")
            for i in range(0,len(context_seq),60):
                outf_context.write(context_seq[i:i+60]+"\n")
            # wt1
            outf_context.write(">" + context_id + "_wt1\n")
            for i in range(0,len(wt1_context_seq),60):
                outf_context.write(wt1_context_seq[i:i+60]+"\n")
            # wt2
            outf_context.write(">" + context_id + "_wt2\n")
            for i in range(0,len(wt2_context_seq),60):
                outf_context.write(wt2_context_seq[i:i+60]+"\n")
        if genome_size < 1000000:
            outf_context.write(">" + "Junk Sequence_" + str(1000000-genome_size) + "bp\n")
            for i in range(0,(1000000-genome_size)/60):
                outf_context.write("N"*60 + "\n")
        outf_context.close()

    def generate_bed(self):
#        outf = open(self.bed,"w")
#        outf_vis = open(self.bed+".visualisation.bed","w")

        context_map = {}
#        print self.fusion_table
        with open(self.fusions_table) as csvfile:
            inf = csv.reader(csvfile, delimiter=';')
            header = inf.next()
#            header = inf.next()
#            for row in inf:
#            next(inf)
            col = {}
            for i,colname in enumerate(header):
                col[colname] = i

            for c,row in enumerate(inf):
#                elements = line.rstrip().split(";")
                ftid = row[col["FTID"]]
                fgid = row[col["FGID"]]
                context_sequence = row[col["context_sequence"]]
                wt1_context_sequence = row[col["wt1_context_sequence"]]
                wt2_context_sequence = row[col["wt2_context_sequence"]]
                bp1_split = row[col["Breakpoint1"]].split(":")
                bp2_split = row[col["Breakpoint2"]].split(":")

#                print fgid
                genes = []
                gene1 = ""
                gene2 = ""
                try:
                    genes = fgid.split("_")
                    gene1 = genes[0]
                    gene2 = genes[2]
                except:
                    genes = fgid.split("_")
                    gene1 = genes[0]
                    gene2 = genes[1]
#                if gene2 == "intronic":
#                    print "INTRONIC!!!"
#                transcript_sequence = row[col["transcript_sequence"]]
                
                context_id = fgid + "_" + row[col["context_sequence_id"]]
                rel_bp_ft = row[col["context_sequence_bp"]]
                rel_bp_wt1 = row[col["wt1_context_sequence_bp"]]
                rel_bp_wt2 = row[col["wt2_context_sequence_bp"]]
#                rel_bp2 = row[col["relative_breakpoint2"]]
#                rel_bp1_ctx = row[col["relative_breakpoint1_ctx"]]
#                rel_bp2_ctx = row[col["relative_breakpoint2_ctx"]]

#                id = ftid + "_fg_" + str(c)
#                id_wt1 = ftid + "_wt1_" + str(c)
#                id_wt2 = ftid + "_wt2_" + str(c)

#                print id
#                print context_sequence
#                print "Rel bp1: " + rel_bp1
#                print "Rel bp2: " + rel_bp2
#                print "Rel bp1 ctx: " + rel_bp1_ctx
#                print "Rel bp2 ctx: " + rel_bp2_ctx
                location_ft = location_wt1 = location_wt2 = 0
                if rel_bp_ft != "NA":
                    location_ft = int(rel_bp_ft)
                if rel_bp_wt1 != "NA":
                    location_wt1 = int(rel_bp_wt1)
                if rel_bp_wt2 != "NA":
                    location_wt2 = int(rel_bp_wt2)
#                elif rel_bp1 == "NA" and rel_bp2 != "NA":
#                    location = int(rel_bp2)
#                elif rel_bp1 != "NA" and rel_bp2 == "NA":
#                    location = int(rel_bp1)
#                elif rel_bp1 == "NA" and rel_bp2 == "NA" and rel_bp1_ctx != "NA":
#                    location = transcript_sequence.find(context_sequence) + int(rel_bp1_ctx)
#                elif rel_bp1 == "NA" and rel_bp2 == "NA" and rel_bp1_ctx == "NA" and rel_bp2_ctx != "NA":
#                    location = transcript_sequence.find(context_sequence) + int(rel_bp2_ctx)
#                elif rel_bp1 == "NA" and rel_bp2 == "NA" and rel_bp1_ctx == "NA" and rel_bp2_ctx == "NA":
#                    location = transcript_sequence.find(context_sequence) + 100


#                chrom1 = id
#                start1 = 0
#                end1 = location_ft
                strand1 = ""
                try:
                    strand1 = bp1_split[2]
                except:
                    strand1 = "+"
#                chrom2 = id
#                start2 = location_ft
#                end2 = len(transcript_sequence)
                strand2 = ""
                try:
                    strand2 = bp2_split[2]
                except:
                    strand2 = "-"
                context_map[context_id] = [location_ft,location_wt1,location_wt2,gene1,strand1,gene2,strand2,len(context_sequence),len(wt1_context_sequence),len(wt2_context_sequence)]
#            print context_map[context_id]
                # ft
#                outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1,start1,end1,gene1,"",strand1,start1,end1,"0",1,end1-start1,start1))
#                outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom2,start2,end2,gene2,"",strand2,start2,end2,"0",1,end2-start2,start2))
                # wt1
#                outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1,start1,end1,gene1,"",strand1,start1,end1,"0",1,end1-start1,start1))
#                outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom2,start2,end2,gene2,"",strand2,start2,end2,"0",1,end2-start2,start2))
                # wt2
#                outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1,start1,end1,gene1,"",strand1,start1,end1,"0",1,end1-start1,start1))
#                outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom2,start2,end2,gene2,"",strand2,start2,end2,"0",1,end2-start2,start2))
                # ft
#                outf_vis.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1,start1,(end1-5),gene1,"",strand1,start1,(end1-5),"0",1,(end1-5)-start1,start1))
#                outf_vis.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom2,(start2+5),end2,gene2,"",strand2,(start2+5),end2,"0",1,end2-(start2+5),(start2+5)))
                # wt1
#                outf_vis.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1,start1,(end1-5),gene1,"",strand1,start1,(end1-5),"0",1,(end1-5)-start1,start1))
#                outf_vis.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom2,(start2+5),end2,gene2,"",strand2,(start2+5),end2,"0",1,end2-(start2+5),(start2+5)))
                # wt2
#                outf_vis.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom1,start1,(end1-5),gene1,"",strand1,start1,(end1-5),"0",1,(end1-5)-start1,start1))
#                outf_vis.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom2,(start2+5),end2,gene2,"",strand2,(start2+5),end2,"0",1,end2-(start2+5),(start2+5)))
#            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom1,start1,end1,chrom2,start2,end2,fusion_gene,"",strand1,strand2,"")
#        outf.close()
#        outf_vis.close()

        outf_context = open(self.bed,"w")
        outf_context_vis = open(self.bed+".visualisation.bed","w")

        for context_id in context_map:
            elements = context_map[context_id]
#        print elements
            start_ft = 0
            start_wt1 = 0
            start_wt2 = 0
            end1_ft = elements[0]
            end1_wt1 = elements[1]
            end1_wt2 = elements[2]
            gene1 = elements[3]
            strand1 = elements[4]
#            start2_ft = end1
#            start2_wt2 = end_

            gene2 = elements[5]
            strand2 = elements[6]
            end2_ft = elements[7]
            end2_wt1 = elements[8]
            end2_wt2 = elements[9]

            # ft
            outf_context.write("%s_ft\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,start_ft,end1_ft,gene1,"",strand1,start_ft,end1_ft,"0",1,end1_ft,start_ft))
            outf_context.write("%s_ft\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,end1_ft,end2_ft,gene2,"",strand2,end1_ft,end2_ft,"0",1,end2_ft-end1_ft,end1_ft))
            # wt1
            outf_context.write("%s_wt1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,start_wt1,end1_wt1,gene1,"",strand1,start_wt1,end1_wt1,"0",1,end1_wt1,start_wt1))
            outf_context.write("%s_wt1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,end1_wt1,end2_wt1,gene2,"",strand2,end1_wt1,end2_wt1,"0",1,end2_wt1-end1_wt1,end1_wt1))
            # wt2
            outf_context.write("%s_wt2\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,start_wt2,end1_wt2,gene1,"",strand1,start_wt2,end1_wt2,"0",1,end1_wt2,start_wt2))
            outf_context.write("%s_wt2\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,end1_wt2,end2_wt2,gene2,"",strand2,end1_wt2,end2_wt2,"0",1,end2_wt2-end1_wt2,end1_wt2))

            # ft
            outf_context_vis.write("%s_ft\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,start_ft,(end1_ft-5),gene1,"",strand1,start_ft,(end1_ft-5),"0",1,(end1_ft-5),start_ft))
            outf_context_vis.write("%s_ft\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,(end1_ft+5),end2_ft,gene2,"",strand2,(end1_ft+5),end2_ft,"0",1,end2_ft-(end1_ft+5),(end1_ft+5)))
            # wt1
            outf_context_vis.write("%s_wt1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,start_wt1,(end1_wt1-5),gene1,"",strand1,start_wt1,(end1_wt1-5),"0",1,(end1_wt1-5),start_wt1))
            outf_context_vis.write("%s_wt1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,(end1_wt1+5),end2_wt1,gene2,"",strand2,(end1_wt1+5),end2_wt1,"0",1,end2_wt1-(end1_wt1+5),(end1_wt1+5)))
            # wt2
            outf_context_vis.write("%s_wt2\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,start_wt2,(end1_wt2-5),gene1,"",strand1,start_wt2,(end1_wt2-5),"0",1,(end1_wt2-5),start_wt2))
            outf_context_vis.write("%s_wt2\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (context_id,(end1_wt2+5),end2_wt2,gene2,"",strand2,(end1_wt2+5),end2_wt2,"0",1,end2_wt2-(end1_wt2+5),(end1_wt2+5)))
        outf_context.close()
        outf_context_vis.close()

    def build_index(self):
        star_idx_path = os.path.join(self.working_dir,"STAR_idx")
#        hisat2_idx_path = os.path.join(self.working_dir,"hisat2_index","fusions").rstrip("/")
#        bowtie_idx_path = os.path.join(self.working_dir,"bowtie_index","fusions").rstrip("/")
        
#        for folder in [star_idx_path,hisat2_idx_path,bowtie_idx_path]:
        for folder in [star_idx_path]:
            create_folder(folder)

        shell_script = os.path.join(self.working_dir,"generate_index.sh")
        out_shell = open(shell_script,"w")
#        fasta_file = self.fasta
        sa_index_nbases = int(math.log(self.calc_genome_size()) / 2 - 1)
        cmd_star = "%s --runMode genomeGenerate --runThreadN 12 --genomeSAindexNbases %s --genomeDir %s --genomeFastaFiles %s" % (self.cfg.get('commands','star_cmd'),sa_index_nbases,star_idx_path,self.fasta)
#        cmd_hisat2 = "/kitty/code/hisat2-2.0.4/hisat2-build -p 12 %s %s" % (fasta_file,hisat2_idx_path)
#        cmd_bowtie = "/kitty/code/bowtie-1.1.2-linux-x86_64/bowtie-build %s %s" % (fasta_file,bowtie_idx_path)

        out_shell.write("#!/bin/sh\n\n")
        out_shell.write("working_dir="+self.working_dir+"\n")
        out_shell.write("echo \"Generating STAR index\"\n")
        out_shell.write(cmd_star+"\n")
#        out_shell.write("echo \"Generating Hisat2 index\"\n")
#        out_shell.write(cmd_hisat2+"\n")
#        out_shell.write("echo \"Generating Bowtie index\"\n")
#        out_shell.write(cmd_bowtie+"\n")
        out_shell.write("echo \"Processing done!\"\n")
        out_shell.close()

#        wd_id = os.path.basename(os.path.normpath(self.working_dir)).upper()
#        index_id = 0
#        for i in range(len(index_base)):
#            index_id += int(ord(self.working_dir[i]))

        if not os.path.exists(os.path.join(star_idx_path,"Genome")):
#            self.submit_job("STAR_IDX",wd_id,cmd_star,12,10,os.path.join(self.working_dir,"STAR_idx"),"")
            self.submit_job_nonqueue(cmd_star,star_idx_path)
#        if not os.path.exists(os.path.join(self.working_dir,"hisat2_index","fusions.1.ht2")):
#            self.submit_job("HISAT2_IDX",wd_id,cmd_hisat2,12,10,os.path.join(self.working_dir,"hisat2_index"),"")
#        if not os.path.exists(os.path.join(self.working_dir,"bowtie_index","fusions.1.ebwt")):
#            self.submit_job("BOWTIE_IDX",wd_id,cmd_bowtie,12,10,os.path.join(self.working_dir,"bowtie_index"),"")

    def run(self):
        '''This function starts processing of the samples.'''
        fastqs = self.get_fastq_files(self.input)
        file_array = sorted(fastqs)
#        print file_array
        left = []
        right = []
        for i in range(len(file_array)):
            file_array_split = file_array[i].rsplit("_",2)
            skewer_split = ""
            try:
                skewer_split = file_array[i].rsplit("-",1)[1].rsplit(".",2)[0]
            except:
                skewer_split = ""
            if "R1" in file_array_split[1] or skewer_split == "pair1":
                left.append(file_array[i])
            elif "R2" in file_array_split[1] or skewer_split == "pair2":
                right.append(file_array[i])
        for i,ele in enumerate(left):
            if len(right) != 0:
                print "Sample 1: " + left[i]
                print "Sample 2: " + right[i]
                self.execute_pipeline(left[i],right[i])

    def execute_pipeline(self,file_1,file_2):

        star_genome_path = os.path.join(self.working_dir,"STAR_idx")

        star_path = os.path.join(self.working_dir,"star")

        for folder in [star_path]:
            IOMethods.create_folder(folder)

        sam_file = os.path.join(star_path,"Aligned.out.sam")
        bam_file = os.path.join(star_path,"Aligned.out.bam")
        
        cmd_sambam_star = "%s view -bS %s | %s sort -o %s && %s index %s && rm %s" % (self.cfg.get('commands','samtools_cmd'),sam_file,self.cfg.get('commands','samtools_cmd'),bam_file,self.cfg.get('commands','samtools_cmd'),bam_file,sam_file)

        cmd_star = "%s --genomeDir %s --readFilesCommand 'gzip -d -c -f' --readFilesIn %s %s --outSAMmode Full --outFilterMultimapNmax 100 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 --runThreadN 12" % (self.cfg.get('commands','star_cmd'),star_genome_path,file_1,file_2)

        cmd_class_add = "%s -i %s -b %s -o %s -m add" % (self.cfg.get('commands','classification_cmd'),bam_file,os.path.join(self.working_dir,self.bed),os.path.join(self.working_dir,"Classification.csv"))

        if not os.path.exists(os.path.join(star_path,"Log.final.out")):
            cpu = 6
            mem = 20
            self.submit_job("STAR", sample_id," && ".join([cmd_star, cmd_sambam_star]), cpu, mem, star_path, dependencies, "")

        if not os.path.exists(os.path.join(self.working_dir,"Classification.csv")):
            cpu = 1
            mem = 80
            self.submit_job("QUANT", sample_id, cmd_class_add, cpu, mem, output_results_folder, Queueing.get_jobs_with_name(sample_id), "")

        print("Processing complete for {}".format(self.working_dir)

    def submit_job(self, module, sample_id, cmd, cores, mem_usage, output_results_folder, dependencies, sched):
        '''This function submits a job with the corresponding resources to slurm.'''
        job_name = "-".join([module,sample_id])
        already_running = Queueing.get_jobs_with_name(job_name)
        if len(already_running) == 0:
            Queueing.submit(job_name, cmd, cores, mem_usage, output_results_folder, dependencies, self.partitions, sched)
            time.sleep(3)
        else:
            print(job_name, " already running!")

        return job_name
       

def main():
    parser = ArgumentParser(description='Processing of demultiplexed FASTQs')

    parser.add_argument('-i', '--input', dest='input', nargs='+', help='Specify the fastq folder(s) or fastq file(s) to process.',required=True)
    parser.add_argument('-f', '--fusions', dest='fusions', help='Specify the reference fusions table to base index on.',required=True)
#    parser.add_argument('-m','--min-read-length-perc',dest='min_read_length_perc',type=float,help='Specify minimum read length percentage',default=0.75)
#    parser.add_argument('--no-qc', dest='no_qc', action='store_true', help='Specify if QC shall be skipped on samples')
    parser.add_argument('-o', '--output-folder', dest='output_folder', help='Specify the folder to save the results into.',required=True)

    parser.add_argument('-p','--partitions',dest='partitions',help='Comma-separated list of partitions to use for queueing.',default='prod')
    args = parser.parse_args()

    cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)),'config.ini'))

    script_call = "python " + os.path.realpath(__file__) + " " + " ".join(sys.argv[1:])
    
    ct = CustomTranscriptome(cfg,args.input,args.fusions,args.output_folder,args.partitions)
#    ct.csv_to_fasta()
#    ct.build_index()
    ct.run()

    outf = open(os.path.join(args.output_folder,"custom_transcriptome.sh"),"w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

if __name__ == '__main__':
    main()
