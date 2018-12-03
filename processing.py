#!/usr/bin/env python

"""
Main class for processing large cohorts of data with easyfuse pipeline
"""

from __future__ import print_function
import os
import sys
import time
from argparse import ArgumentParser

from misc.config import Config
from misc.samples import Samples
import misc.queue as Queueing
import misc.io_methods as IOMethods

class Processing(object):
    def __init__(self, config, input_paths, working_dir, ref_genome, min_rl_perc, filter_reads, partitions):
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
        self.filter_reads = filter_reads
        self.partitions = partitions
        self.sched = self.cfg.get('general', 'queueing_system')

    def run(self):
        '''This function starts processing of the samples.'''

        if icam_run:
            self.logger.info("Started processing of icam data; searching fastq files in \"{}\"".format(self.input_path))
            # print(len(self.samples.sample_map)) <- 0 for empty list
            _, fastq_dict, counter_list = urlaio.get_icam_reads_data(self.input_path, list(self.samples.sample_map), self.logger)
            count_all = sum(counter_list)
            self.logger.debug("Found a total of {0} datasets in the specified folder, {1} have been processed previously, {2} are new.".format(count_all, counter_list[0], counter_list[1]))
            self.logger.debug("{0} datasets could not be processed due to an unexpected format of the demux_stats.csv file.".format(counter_list[2]))
            self.logger.debug("{0} datasets were not be processed because their read count was below 75M.".format(counter_list[3]))
            self.logger.debug("{0} datasets could not be processed due problems in the fastq file localization.".format(counter_list[4]))
            for i, sample_id_key in enumerate(sorted(fastq_dict)):
                sample_id = "_".join(sample_id_key.split("_")[:2])
                self.logger.info("Processing Sample ID: {} (paired end)".format(sample_id))
                self.logger.info("Sample 1: {}".format(fastq_dict[sample_id_key][0]))
                self.logger.info("Sample 2: {}".format(fastq_dict[sample_id_key][1]))
                self.execute_pipeline(fastq_dict[sample_id_key][0], fastq_dict[sample_id_key][1], sample_id, ref_genome, ref_trans, tool_num_cutoff, filter_reads)
        else:
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

    def execute_pipeline(self, fq1, fq2, sample_id):
        '''This function executes a defined set of tools on a given sample and creates the respective folder structure.'''
        
        organism = "human" if self.ref_genome in ['hg18', 'hg19', 'hg38'] else "mouse"

        self.samples.add_sample(sample_id, "NA", "R1:{};R2:{}".format(fq1, fq2))
        star_genome_path = self.cfg.get('star_indexes', 'star_index_{}'.format(self.ref_genome))
        starfusion_genome_path = self.cfg.get('star_fusion_indexes', 'star_fusion_index_{}'.format(self.ref_genome))
        fusioncatcher_genome_path = self.cfg.get('fusioncatcher_indexes', 'fusioncatcher_index_{}'.format(self.ref_genome))

        # Path' to specific indices
        bowtie_index_path = self.cfg.get("indices", "{0}_bowtie_{1}".format(ref_trans, ref_genome))
        star_index_path = self.cfg.get("indices", "{0}_star_{1}_sjdb{2}".format(ref_trans, ref_genome, 49))
        #        star_index_path = self.cfg.get("indices", "{0}_star_{1}_trans".format(ref_trans, ref_genome))
        kallisto_index_path = self.cfg.get("indices", "{0}_kallisto_{1}".format(ref_trans, ref_genome))
        pizzly_cache_path = "{}.pizzlyCache.txt".format(genes_gtf_path)
        starfusion_index_path = self.cfg.get("indices", "{0}_starfusion_{1}".format(ref_trans, ref_genome))
        infusion_cfg_path = self.cfg.get("otherFiles", "{0}_infusion_cfg_{1}".format(ref_trans, ref_genome))
        starchip_param_path = self.cfg.get("otherFiles", "{0}_starchip_param_{1}".format(ref_trans, ref_genome))


#        output_results_folder = os.path.join(self.working_dir, "Sample_{}".format(sample_id), "scratch")
        output_results_folder = os.path.join(self.working_dir, sample_id, "scratch")
        filtered_reads_path = os.path.join(output_results_path, "filtered_reads")
        expression_path = os.path.join(output_results_folder, "expression")
        kallisto_path = os.path.join(expression_path, "kallisto")
        star_path = os.path.join(expression_path, "star")
                
        fastqc_path = os.path.join(output_results_folder, "fastqc")
        skewer_path = os.path.join(output_results_folder, "skewer")
        qc_table_path = os.path.join(output_results_folder, "qc_table.txt")
        overrepresented_path = os.path.join(output_results_folder, "overrepresented.fa")

        fusion_path = os.path.join(output_results_folder, "fusion")
        starfusion_path = os.path.join(fusion_path, "starfusion")
        mapsplice_path = os.path.join(fusion_path, "mapsplice")
        fusioncatcher_path = os.path.join(fusion_path, "fusioncatcher")
        soapfuse_path = os.path.join(fusion_path, "soapfuse")
        infusion_path = os.path.join(fusion_path, "infusion")

        # IVAC specific tools
        pizzly_path = os.path.join(fusion_path, "pizzly")
        starchip_path = os.path.join(fusion_path, "starchip")

        fetchdata_path = os.path.join(self.working_dir, sample_id, "fetchdata")

        trimmed = os.path.join(skewer_path, "out_file-trimmed.fastq.gz")
        trimmed_1 = os.path.join(skewer_path, "out_file-trimmed-pair1.fastq.gz")
        trimmed_2 = os.path.join(skewer_path, "out_file-trimmed-pair2.fastq.gz")

        unmapped_1 = os.path.join(output_results_folder, "Unmapped.out.mate1")
        unmapped_2 = os.path.join(output_results_folder, "Unmapped.out.mate2")

        potential_1 = os.path.join(filter_path, "fastqs", "potential_R1_001.fastq.gz")
        potential_2 = os.path.join(filter_path, "fastqs", "potential_R2_001.fastq.gz")


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
                       infusion_path,
                       pizzly_path,
                       starchip_path,
                       fetchdata_path
            ]:
            IOMethods.create_folder(folder)


        cmd_fastqc = "{} --nogroup --extract -t 6 -o {} {} {}".format(self.cfg.get('commands', 'fastqc_cmd'), fastqc_path, fq1, fq2)
        cmd_readcoon = "{} -i {}/*/fastqc_data.txt -q {} -o {}".format(self.cfg.get('commands', 'readcoon_cmd'), fastqc_path, qc_table_path, overrepresented_path)
        cmd_skewer = "{} -q {} -m {} -i {} {} -o {}".format(self.cfg.get('commands', 'skewer_wrapper_cmd'), qc_table_path, self.min_rl_perc, fq1, fq2, skewer_path)


        # (0) Readfilter
        cmd_star_filter = "{0} --genomeDir {1} --outFileNamePrefix {2} --readFilesCommand zcat --readFilesIn {3} {4} --outFilterMultimapNmax 100 --outSAMmultNmax 1 --chimSegmentMin 10 --chimJunctionOverhangMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --seedSearchStartLmax 20 --winAnchorMultimapNmax 50 --outSAMtype BAM Unsorted --chimOutType Junctions WithinBAM --outSAMunmapped Within KeepPairs --runThreadN waiting_for_cpu_number".format(self.cfg.get('commands', 'star_cmd'), star_index_path, os.path.join(filtered_reads_path, "{}_".format(sample_id)), trimmed_1, trimmed_2)

        cmd_read_filter = "{0} --input {1}_Aligned.out.bam --output {1}_Aligned.out.filtered.bam".format(self.cfg.get('commands', 'readfilter_cmd'), os.path.join(filtered_reads_path, sample_id))

        # re-define fastq's if filtering is on
        fq0 = ""
        if self.filter_reads:
            fq0 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(".fastq.gz", "_filtered_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered.fastq.gz"))
            tools.insert(0, "Readfilter")

        cmd_bam_to_fastq = "{0} fastq -0 {1} -1 {2} -2 {3} {4}_Aligned.out.filtered.bam".format(self.cfg.get('commands', 'samtools_cmd'), fq0, fq1, fq2, os.path.join(filtered_reads_path, sample_id))


#        cmd_filter = "{} -i {} {} -o {} -g {}".format(self.cfg.get('commands', 'filter_cmd'), trimmed_1, trimmed_2, filter_path, self.ref_genome)

#        if self.filter_reads:
#            trimmed_1 = potential_1
#            trimmed_2 = potential_2        

        cmd_starfusion = "{} --genome_lib_dir {} --left_fq {} --right_fq {} --output_dir {}".format(self.cfg.get('commands', 'starfusion_cmd'), starfusion_genome_path, trimmed_1, trimmed_2, starfusion_path)

        cpu, mem = self.cfg.get('resources','fusioncatcher').split(",")
        cmd_fusioncatcher = "{} -d {} --input {},{} --output {} -p {}".format(self.cfg.get('commands', 'fusioncatcher_cmd'), fusioncatcher_genome_path, trimmed_1, trimmed_2, fusioncatcher_path, cpu)

        cmd_soapfuse = "{} -c {} -q {} -i {} -o {} -g {}".format(self.cfg.get('commands', 'soapfuse_wrapper_cmd'), self.cfg._file, qc_table_path, " ".join([trimmed_1, trimmed_2]), soapfuse_path, organism)

        cpu, mem = self.cfg.get('resources', 'mapsplice').split(",")
        cmd_mapsplice = "{} -c {} -i {} -t {} -o {} -g {}".format(self.cfg.get('commands','mapsplice_wrapper_cmd'), self.cfg._file, " ".join([trimmed_1, trimmed_2]), cpu, mapsplice_path, self.ref_genome)

        cmd_infusion = "{} -1 {} -2 {} --skip-finished --min-unique-alignment-rate 0 --min-unique-split-reads 0 --allow-non-coding --out-dir {} {}".format(self.cfg.get('commands', 'infusion_cmd'), trimmed_1, trimmed_2, infusion_path, self.cfg.get('commands', 'infusion_cfg_{}'.format(organism)))

        cmd_fetchdata = "{} -i {}".format(self.cfg.get('commands', 'fetchdata_cmd'), output_results_folder)

        # (10) De novo assembly of fusion transcripts
        # urla: This is currently still under active development and has not been tested thoroughly
        cmd_denovoassembly = "{0} -i waiting_for_gene_list_input -c {1} -b {2}_Aligned.out.bam -g {3} -t {4} -o waiting_for_assembly_out_dir".format(self.cfg.get('commands', 'denovoassembly_cmd'), self.cfg.get_path(), os.path.join(filtered_reads_path, sample_id), ref_genome, ref_trans)

        cmd_samples = "{} --output-file={} --sample-id={} --action=append_state --state='".format(self.cfg.get('commands', 'samples_cmd'), self.samples.infile, sample_id)

        state_tools = self.samples.get_tool_list_from_state(sample_id)
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

            if not os.path.exists(potential_1) and self.filter_reads:
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


            if not self.filter_reads:
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


        exe_tools = [
            "FastQC",
            "Readcoon",
            "Skewer",
            "ReadFilter",
            "Kallisto",
            "Star",
            "Mapsplice",
            "Fusioncatcher",
            "Starfusion",
            "Soapfuse",
            "Infusion",
            "Pizzly",
            "Starchip",
            "Fetchdata",
            "Assembly"
        ]

        exe_cmds = [
            cmd_fastqc,
            cmd_readcoon,
            cmd_skewer,
            " && ".join([cmd_star_filter, cmd_read_filter, cmd_bam_to_fastq]), #0
            cmd_kallisto, #1
            cmd_star, #2
            " && ".join([cmd_extr_fastq1, cmd_extr_fastq2, cmd_mapsplice]), #3
            cmd_fusioncatcher, #4
            cmd_starfusion, #5
            " && ".join([cmd_pizzly, cmd_pizzly2]), #6
            cmd_starchip, #7
            cmd_infusion, #8
            cmd_fetchdata, #9
            cmd_denovoassembly #10
        ]

        exe_path = [
            fastqc_path,
            output_results_folder,
            skewer_path,
            filtered_reads_path, #0
            kallisto_path, #1
            star_path, #2
            mapsplice_path, #3
            fusioncatcher_path, #4
            starfusion_path, #5
            pizzly_path, #6
            starchip_path, #7
            infusion_path, #8
            fetchdata_path, #9
            "" #10
        ]

        # create and submit slurm job if the tool is requested and hasn't been run before
        for i, tool in enumerate(exe_tools, 0):
            if tool in tools:
                dependency = []
                if tool == "Pizzly" and "Kallisto" not in tools and "Kallisto" not in state_tools:
                    self.logger.error("Pizzly was selected w/o Kallisto and w/o existing Kallisto data -> this does not work!")
                    sys.exit(99)
                elif (tool == "Starfusion" or tool == "Starchip") and "Star" not in tools and "Star" not in state_tools:
                    self.logger.error("Starchip and/or Starfusion was selected w/o Star and w/o existing Star data -> this does not work!")
                    sys.exit(99)
                if tool in state_tools:
                    # urla: the primary idea behind this flag is to allow multiple fetchdata executions during processing
                    #       nevertheless, re-processing of the same data with a newer version of a tool will also be straightforward (but overwriting previous results, of course)
                    if self.overwrite:
                        self.logger.info("Executing {0} although it looks like a previous run finished successfully. Results in {1} may be overwritten".format(tool, exe_path[i]))
                    else:
                        self.logger.info("Skipping {0} as it looks like a previous run finished successfully. Results should be in {1}".format(tool, exe_path[i]))
                        continue
                self.logger.debug("Submitting {} run to slurm".format(tool))
                cpu, mem = self.cfg.get("resources", tool.lower()).split(",")
                uid = "-".join([tool, sample_id, str(int(round(time.time())))])
                if tool == "Star":
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string", "{}/{}_".format(exe_path[i], sample_id)).replace("waiting_for_cpu_number", str(cpu))
                else:
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string", exe_path[i]).replace("waiting_for_cpu_number", str(cpu))
                cmd = " && ".join([exe_cmds[i], cmd_samples + uid + "'"])
                # Managing dependencies
                if tool == "Pizzly":
                    dependency = Queueing.get_jobs_with_name_slurm("Kallisto-{}".format(sample_id))
                elif tool == "Starfusion" or tool == "Starchip":
                    dependency = Queueing.get_jobs_with_name_slurm("Star-{}".format(sample_id))
                elif tool == "Fetchdata":
                    dependency = Queueing.get_jobs_with_name_slurm(sample_id)
                elif tool == "Assembly":
                    dependency = Queueint.get_jobs_with_name_slurm("Fetchdata-{}".format(sample_id))
                if self.filter_reads:
                    dependency.extend(Queueing.get_jobs_with_name_slurm("Readfilter-{}".format(sample_id)))
                self.logger.debug("Submitting slurm job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd, exe_path[i], dependency))
                self.submit_job(uid, cmd, cpu, mem, exe_path[i], dependency)
            else:
                self.logger.info("Skipping {0} as it is not selected for execution (Selected are: {1})".format(tool, tools))
                                    

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
    parser.add_argument('--filter-reads', dest='filter_reads', action='store_true', help='Select this option if you do want to pre-filter your input reads.')
    parser.add_argument('-p','--partitions', dest='partitions', help='Comma-separated list of partitions to use for queueing.', default='prod')
    args = parser.parse_args()

    if args.config == "":
        cfg = Config(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini'))
    else:
        cfg = Config(args.config)

    script_call = "python {} {}".format(os.path.realpath(__file__), " ".join(sys.argv[1:]))

    proc = Processing(cfg, args.input, args.output_folder, args.ref_genome, args.min_read_length_perc, args.filter_reads, args.partitions)
    proc.run()

    outf = open(os.path.join(args.output_folder, "process.sh"),"w")
    outf.write("#!/bin/sh\n\n")
    outf.write(script_call)
    outf.close()

if __name__ == '__main__':
    main()
