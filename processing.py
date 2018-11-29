#!/usr/bin/env python

"""
Main class for processing large cohorts of data with easyfuse pipeline
"""

import os
import sys
import time
from argparse import ArgumentParser

from misc.config import Config
from misc.samples import Samples
import misc.queue as Queueing
import misc.io_methods as IOMethods

class Processing(object):
    def __init__(self, config, input_paths, working_dir, ref_genome, min_rl_perc, no_filter, partitions):
        '''This function inits the input parameters and creates the working directory.'''
        self.cfg = config
        self.input_paths = []
        for i in input_paths:
            self.input_paths.append(i.rstrip("/"))
        self.working_dir = working_dir.rstrip("/")
        IOMethods.create_folder(self.working_dir)
        self.samples = Samples(os.path.join(self.working_dir, "samples.csv"))
        self.ref_genome = ref_genome
        self.min_rl_perc = min_rl_perc
        self.no_filter = no_filter
        self.partitions = partitions
        self.sched = self.cfg.get('general', 'queueing_system')

    def run(self):
        '''This function starts processing of the samples.'''
        fastqs = IOMethods.get_fastq_files(self.input_paths)
        file_array = sorted(fastqs)
        left = []
        right = []
        for i in range(len(file_array)):
            file_array_split = file_array[i].rsplit("_",2)
            if "R1" in file_array_split[1]:
                left.append(file_array[i])
            elif "R2" in file_array_split[1]:
                right.append(file_array[i])
        for i,ele in enumerate(left):
            sample_id = left[i].rsplit("/",1)[1].rsplit("_",2)[0]
            if len(right) != 0:
                print("Processing Sample ID: {} (paired-end)".format(sample_id))
                print("Sample 1: {}".format(left[i]))
                print("Sample 2: {}".format(right[i]))
                self.execute_pipeline(left[i], right[i], sample_id)

    def execute_pipeline(self, file_1, file_2, sample_id):
        '''This function executes a defined set of tools on a given sample and creates the respective folder structure.'''
        
        organism = "human" if self.ref_genome in ['hg18', 'hg19', 'hg38'] else "mouse"

        self.samples.add(sample_id,"NEW [0/4]")
        star_genome_path = self.cfg.get('star_indexes', 'star_index_{}'.format(self.ref_genome))
        starfusion_genome_path = self.cfg.get('star_fusion_indexes', 'star_fusion_index_{}'.format(self.ref_genome))
        fusioncatcher_genome_path = self.cfg.get('fusioncatcher_indexes', 'fusioncatcher_index_{}'.format(self.ref_genome))

        output_results_folder = os.path.join(self.working_dir, "Sample_{}".format(sample_id), "scratch")
        expression_path = os.path.join(output_results_folder, "expression")
        fastqc_path = os.path.join(output_results_folder, "fastqc")
        skewer_path = os.path.join(output_results_folder, "skewer")
        filter_path = os.path.join(output_results_folder, "filtered")
        hla_path = os.path.join(output_results_folder, "seq2hla")
        isoform_path = os.path.join(output_results_folder, "isoform")
        fusion_path = os.path.join(output_results_folder, "fusion")
        starfusion_path = os.path.join(fusion_path, "starfusion")
        mapsplice_path = os.path.join(fusion_path, "mapsplice")
        fusioncatcher_path = os.path.join(fusion_path, "fusioncatcher")
        soapfuse_path = os.path.join(fusion_path, "soapfuse")
        infusion_path = os.path.join(fusion_path, "infusion")
        qc_table_path = os.path.join(output_results_folder, "qc_table.txt")
        overrepresented_path = os.path.join(output_results_folder, "overrepresented.fa")

        for folder in [output_results_folder, 
                       expression_path, 
                       fastqc_path, 
                       skewer_path, 
                       filter_path, 
                       hla_path, 
                       fusion_path, 
                       starfusion_path, 
                       fusioncatcher_path, 
                       soapfuse_path, 
                       mapsplice_path, 
                       infusion_path
            ]:
            IOMethods.create_folder(folder)

        trimmed = os.path.join(skewer_path, "out_file-trimmed.fastq.gz")
        trimmed_1 = os.path.join(skewer_path, "out_file-trimmed-pair1.fastq.gz")
        trimmed_2 = os.path.join(skewer_path, "out_file-trimmed-pair2.fastq.gz")

        unmapped_1 = os.path.join(output_results_folder, "Unmapped.out.mate1")
        unmapped_2 = os.path.join(output_results_folder, "Unmapped.out.mate2")

        potential_1 = os.path.join(filter_path, "fastqs", "potential_R1_001.fastq.gz")
        potential_2 = os.path.join(filter_path, "fastqs", "potential_R2_001.fastq.gz")

        cmd_samples = "{} --output-file={} --sample-id={} --action=update --state='".format(self.cfg.get('commands', 'samples_cmd'), self.samples.infile, sample_id)

        cmd_fastqc = "{} --nogroup --extract -t 6 -o {} {} {}".format(self.cfg.get('commands', 'fastqc_cmd'), fastqc_path, file_1, file_2)
        cmd_readcoon = "{} -i {}/*/fastqc_data.txt -q {} -o {}".format(self.cfg.get('commands', 'readcoon_cmd'), fastqc_path, qc_table_path, overrepresented_path)
        cmd_skewer = "{} -q {} -m {} -i {} {} -o {}".format(self.cfg.get('commands', 'skewer_wrapper_cmd'), qc_table_path, self.min_rl_perc, file_1, file_2, skewer_path)

        cmd_filter = "{} -i {} {} -o {} -g {}".format(self.cfg.get('commands', 'filter_cmd'), trimmed_1, trimmed_2, filter_path, self.ref_genome)

        if not self.no_filter:
            trimmed_1 = potential_1
            trimmed_2 = potential_2
        cmd_starfusion = "{} --genome_lib_dir {} --left_fq {} --right_fq {} --output_dir {}".format(self.cfg.get('commands', 'starfusion_cmd'), starfusion_genome_path, trimmed_1, trimmed_2, starfusion_path)
        cpu, mem = self.cfg.get('resources','fusioncatcher').split(",")
        cmd_fusioncatcher = "{} -d {} --input {},{} --output {} -p {}".format(self.cfg.get('commands', 'fusioncatcher_cmd'), fusioncatcher_genome_path, file_1, file_2, fusioncatcher_path, cpu)
        cmd_soapfuse = "{} -c {} -q {} -i {} -o {} -g {}".format(self.cfg.get('commands', 'soapfuse_wrapper_cmd'), self.cfg._file, qc_table_path, " ".join([trimmed_1, trimmed_2]), soapfuse_path, organism)
        cpu, mem = self.cfg.get('resources', 'mapsplice').split(",")
        cmd_mapsplice = "{} -c {} -i {} -t {} -o {} -g {}".format(self.cfg.get('commands','mapsplice_wrapper_cmd'), self.cfg._file, " ".join([trimmed_1, trimmed_2]), cpu, mapsplice_path, self.ref_genome)
        cmd_infusion = "{} -1 {} -2 {} --skip-finished --min-unique-alignment-rate 0 --min-unique-split-reads 0 --allow-non-coding --out-dir {} {}".format(self.cfg.get('commands', 'infusion_cmd'), trimmed_1, trimmed_2, infusion_path, self.cfg.get('commands', 'infusion_cfg_{}'.format(organism)))

        cmd_fetchdata = "{} -i {}".format(self.cfg.get('commands', 'fetchdata_cmd'), output_results_folder)

        state = self.samples.get(sample_id)
        if state == "NEW [0/4]":
            if not self.no_qc:
                if not os.path.exists(os.path.join(output_results_folder, "{}_R1_001_fastqc.zip".format(sample_id))):
                    cmd = " && ".join([cmd_fastqc, "{}QC [1/4]'".format(cmd_samples)])
                    cpu,mem = self.cfg.get('resources', 'fastqc').split(",")
                    self.submit_job("FastQC", sample_id, cmd, cpu, mem, output_results_folder, Queueing.get_jobs_by_name(sample_id, self.sched))
        if state in ["QC [1/4]","NEW [0/4]"]:
            if not self.no_qc:
                if not os.path.exists(os.path.join(output_results_folder, "qc_table.txt")):
                    cmd = " && ".join([cmd_readcoon, "{}COONED [2/4]'".format(cmd_samples)])
                    cpu,mem = self.cfg.get('resources', 'readcoon').split(",")
                    self.submit_job("Readcoon", sample_id, cmd, cpu, mem, output_results_folder, Queueing.get_jobs_by_name(sample_id, self.sched))
        if state in ["COONED [2/4]", "QC [1/4]", "NEW [0/4]"]:
            if not self.no_qc:
                if not os.path.exists(os.path.join(skewer_path, "out_file-trimmed-pair1.fastq.gz")):
                    cmd = " && ".join([cmd_skewer, "{}SKEWED [3/4]'".format(cmd_samples)])
                    cpu,mem = self.cfg.get('resources', 'skewer').split(",")
                    self.submit_job("Skewer", sample_id, cmd, cpu, mem, skewer_path, Queueing.get_jobs_by_name(sample_id, self.sched))
        if state in ["SKEWED [3/4]", "COONED [2/4]", "QC [1/4]", "NEW [0/4]"]:

            if not os.path.exists(potential_1) and not self.no_filter:
                cpu,mem = self.cfg.get('resources', 'filter').split(",")
                self.submit_job("FILTER", sample_id, cmd_filter, cpu, mem, filter_path, Queueing.get_jobs_by_name(sample_id, self.sched))

            tools = self.cfg.get('general','tools').split(",")
            dependencies = Queueing.get_jobs_by_name(sample_id)
            if not os.path.exists(os.path.join(starfusion_path, "star-fusion.fusion_candidates.final")) and "starfusion" in tools:
                cpu,mem = self.cfg.get('resources', 'starfusion').split(",")
                self.submit_job("Starfusion", sample_id, cmd_starfusion, cpu, mem, starfusion_path, dependencies)
            if not os.path.exists(os.path.join(mapsplice_path, "fusions_well_annotated.txt")) and "mapsplice" in tools:
                cpu,mem = self.cfg.get('resources','mapsplice').split(",")
                self.submit_job("MapSplice", sample_id, cmd_mapsplice, cpu, mem, mapsplice_path, dependencies)
            if not os.path.exists(os.path.join(fusioncatcher_path, "summary_candidate_fusions.txt")) and "fusioncatcher" in tools:
                cpu,mem = self.cfg.get('resources', 'fusioncatcher').split(",")
                self.submit_job("Fusioncatcher", sample_id, cmd_fusioncatcher, cpu, mem, fusioncatcher_path, dependencies)


            if self.no_filter:
                if not os.path.exists(os.path.join(soapfuse_path, "final_fusion_genes", "Sample_{}".format(sample_id), "Sample_{}.final.Fusion.specific.for.genes".format(sample_id))) and not os.path.exists(os.path.join(soapfuse_path, "final_fusion_genes", "scratch", "scratch.final.Fusion.specific.for.genes")) and "soapfuse" in tools:
                    cpu, mem = self.cfg.get('resources', 'soapfuse').split(",")
                    self.submit_job("Soapfuse", sample_id, cmd_soapfuse, cpu, mem, soapfuse_path, dependencies)
            else:
                if not os.path.exists(os.path.join(soapfuse_path, "final_fusion_genes", "filtered", "filtered.final.Fusion.specific.for.genes")) and not os.path.exists(os.path.join(soapfuse_path, "final_fusion_genes", "scratch", "scratch.final.Fusion.specific.for.genes")) and "soapfuse" in tools:
                    cpu, mem = self.cfg.get('resources', 'soapfuse').split(",")
                    self.submit_job("Soapfuse", sample_id, cmd_soapfuse, cpu, mem, soapfuse_path, dependencies)
            if not os.path.exists(os.path.join(infusion_path, "fusions.detailed.txt")) and "infusion" in tools:
                cpu, mem = self.cfg.get('resources', 'infusion').split(",")
                self.submit_job("Infusion", sample_id, cmd_infusion, cpu, mem, infusion_path, dependencies)

            if not os.path.exists(os.path.join(output_results_folder, "fetchdata_1tool", "Classification.csv")) or not os.path.exists(os.path.join(output_results_folder, "fetchdata_2tool", "Classification.csv")):
                cpu, mem = self.cfg.get('resources', 'fetchdata').split(",")
                self.submit_job("fetchdata", sample_id, cmd_fetchdata, cpu, mem, output_results_folder, Queueing.get_jobs_by_name(sample_id, self.sched))
            cpu, mem = self.cfg.get('resources', 'collect').split(",")
            self.submit_job("Collect", sample_id, "{}DONE [4/4]'".format(cmd_samples), cpu, mem, output_results_folder, Queueing.get_jobs_by_name(sample_id, self.sched))
        if state == "DONE [4/4]":
            print("Processing complete for {}".format(sample_id))

    def submit_job(self, module, sample_id, cmd, cores, mem_usage, output_results_folder, dependencies):
        '''This function submits a job with the corresponding resources to slurm.'''
        job_name = "-".join([module, sample_id])
        already_running = Queueing.get_jobs_by_name(job_name)
        if len(already_running) == 0:
            Queueing.submit(job_name, cmd, cores, mem_usage, output_results_folder, dependencies, self.partitions, self.sched)
            time.sleep(3)
        else:
            print("{} already running!".format(job_name))
        
def main():
    parser = ArgumentParser(description='Fusion Detection on RNA sequencing data')

    parser.add_argument('-i', '--input', dest='input', nargs='+', help='Specify the fastq folder(s) or fastq file(s) to process.', required=True)
    parser.add_argument('-o', '--output-folder', dest='output_folder', help='Specify the folder to save the results into.', required=True)
    parser.add_argument('-c', '--config', dest='config', help='Specify config file.', default="")
    parser.add_argument('-g', '--reference-genome', dest='ref_genome', choices=['hg18','hg19','hg38','mm9','mm10'], help='Specify the reference genome name.', required=True)
    parser.add_argument('-m','--min-read-length-perc', dest='min_read_length_perc', type=float, help='Specify minimum read length percentage', default=0.75)    
    parser.add_argument('--no-filter', dest='no_filter', action='store_true', help='Select this option if you do not want to pre-filter your input reads.')
    parser.add_argument('-p','--partitions', dest='partitions', help='Comma-separated list of partitions to use for queueing.', default='prod')
    args = parser.parse_args()

    if args.config == "":
        cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini'))
    else:
        cfg = Config(args.config)

    script_call = "python {} {}".format(os.path.realpath(__file__), " ".join(sys.argv[1:]))

    proc = Processing(cfg, args.input, args.output_folder, args.ref_genome, args.min_read_length_perc, args.no_filter, args.partitions)
    proc.run()

    outf = open(os.path.join(args.output_folder, "process.sh"),"w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

if __name__ == '__main__':
    main()
