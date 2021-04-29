#!/usr/bin/env python

"""
Read fastq files
create per-sample (sub)folder tree
run gene fusion predictions
start data collection
during all steps, perform slurm scheduling and process monitoring

@author: Tron (PASO), BNT (URLA)
@version: 20190118
"""

import os
import os.path
import sys
import time
import argparse
from shutil import copy
import misc.queueing as Queueing
from misc.samples import SamplesDB
from misc.logger import Logger
import misc.io_methods as IOMethods
from misc.versioncontrol import VersCont
import config as cfg

# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class Processing(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""
    def __init__(self, cmd, input_paths, working_dir):
        """Parameter initiation and work folder creation. Start of progress logging."""
        self.working_dir = os.path.abspath(working_dir)
        self.logger = Logger(os.path.join(self.working_dir, "easyfuse_processing.log"))
        IOMethods.create_folder(self.working_dir, self.logger)
        copy(os.path.join(cfg.module_dir, "config.py"), working_dir)
        self.logger.info("Starting easyfuse: CMD - {}".format(cmd))
        self.input_paths = [os.path.abspath(file) for file in input_paths]
        self.samples = SamplesDB(os.path.join(self.working_dir, "samples.db"))

    # The run method simply greps and organises fastq input files.
    # Fastq pairs (single end input is currently not supported) are then send to "execute_pipeline"
    def run(self, tool_num_cutoff):
        """General parameter setting, identification of fastq files and initiation of processing"""
        self.logger.info("Pipeline Version: {}".format(cfg.version))
        # Checking dependencies
        #VersCont(os.path.join(cfg.module_dir, "dependency_versions.txt")).get_and_print_tool_versions()
        #self.cfg.run_self_test()
        # urla: organism is currently not used for anything, however, this might change; is mouse processing relevant at some point?
        ref_genome = cfg.ref_genome_build
        ref_trans = cfg.ref_trans_version

        self.logger.info("Reference Genome: {0}, Reference Transcriptome: {1}".format(ref_genome, ref_trans))
 #       if self.overwrite:
 #           self.logger.info("#############################################################################")
 #           self.logger.info("")
 #           self.logger.info("Overwrite flag is set => all previously existing results may be overwritten!")
 #           self.logger.info("")
 #           self.logger.info("#############################################################################")


        sample_list = []
        # get fastq files
        left, right, sample_id = IOMethods.get_fastq_files(self.input_paths, self.logger)
        sample_list = sample_id
        for i, _ in enumerate(left):
            if len(left) == len(right):
                self.logger.info("Processing Sample ID: {} (paired end)".format(sample_id[i]))
                self.logger.info("Sample 1: {}".format(left[i]))
                self.logger.info("Sample 2: {}".format(right[i]))
                self.execute_pipeline(left[i], right[i], sample_id[i], ref_genome, ref_trans, tool_num_cutoff)

        
        # summarize all data if selected
        if "Summary" in cfg.tools:
            #dependency = [Queueing.get_jobs_by_name("Fetchdata-{}".format(sample)) for sample in sample_list]
            # urla - note: would be happy to get the dependencies with a stacked LC, but is atm to complicated for me ^^
            dependency = []
            for sample in sample_list:
                dependency.extend(Queueing.get_jobs_by_name("Fetchdata-{}".format(sample), cfg.queueing_system))
            modelling_string = ""
            if cfg.other_files["easyfuse_model"]:
                modelling_string = " --model_predictions"
            cmd_summarize = "python {0} --input {1}{2}".format(os.path.join(cfg.module_dir, "summarize_data.py"), self.working_dir, modelling_string)
            self.logger.debug("Submitting slurm job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd_summarize, self.working_dir, dependency))
            cpu = cfg.resources["summary"]["cpu"]
            mem = cfg.resources["summary"]["mem"]
            self.submit_job("-".join([cfg.pipeline_name, "Summary", str(int(round(time.time())))]), cmd_summarize, cpu, mem, self.working_dir, dependency, cfg.receiver)

    # Per sample, define input parameters and execution commands, create a folder tree and submit runs to slurm
    def execute_pipeline(self, fq1, fq2, sample_id, ref_genome, ref_trans, tool_num_cutoff):
        """Create sample specific subfolder structure and run tools on fastq files"""
        self.samples.add_sample(sample_id, "NA", fq1, fq2)

        refs = cfg.references

        # Genome/Gene references to use
        genome_sizes_path = refs["genome_sizes"]
        genome_chrs_path = refs["genome_fastadir"]
        genes_fasta_path = refs["genes_fasta"]
        genes_gtf_path = refs["genes_gtf"]

        # Path' to specific indices
        indices = cfg.indices
        other_files = cfg.other_files

        bowtie_index_path = indices["bowtie"]
        star_index_path = indices["star"]
#        kallisto_index_path = indices["kallisto"]
#        pizzly_cache_path = "{}.pizzlyCache.txt".format(genes_gtf_path)
        starfusion_index_path = indices["starfusion"]
        fusioncatcher_index_path = indices["fusioncatcher"]
        infusion_cfg_path = other_files["infusion_cfg"]
#        starchip_param_path = other_files["starchip_param"]

        # Output results folder creation - currently included:
        # 1) Gene/Isoform expression: kallisto, star
        # 2) Fusion prediction: mapsplice, pizzly, fusioncatcher, star-fusion, starchip, infusion
        output_results_path = os.path.join(self.working_dir, "Sample_{}".format(sample_id))
        qc_path = os.path.join(output_results_path, "qc")
        skewer_path = os.path.join(qc_path, "skewer")
        qc_table_path = os.path.join(qc_path, "qc_table.txt")
        overrepresented_path = os.path.join(qc_path, "overrepresented.fa")
        filtered_reads_path = os.path.join(output_results_path, "filtered_reads")
        expression_path = os.path.join(output_results_path, "expression")
#        kallisto_path = os.path.join(expression_path, "kallisto")
        star_path = os.path.join(expression_path, "star")
        fusion_path = os.path.join(output_results_path, "fusion")
        mapsplice_path = os.path.join(fusion_path, "mapsplice")
        pizzly_path = os.path.join(fusion_path, "pizzly")
        fusioncatcher_path = os.path.join(fusion_path, "fusioncatcher")
        starfusion_path = os.path.join(fusion_path, "starfusion")
        starchip_path = os.path.join(fusion_path, "starchip")
        infusion_path = os.path.join(fusion_path, "infusion")
        soapfuse_path = os.path.join(fusion_path, "soapfuse")
        fetchdata_path = os.path.join(self.working_dir, "Sample_{}".format(sample_id), "fetchdata")
        fastqc_1 = os.path.join(qc_path, os.path.basename(fq1).rstrip(".fastq.gz") + "_fastqc", "fastqc_data.txt")
        fastqc_2 = os.path.join(qc_path, os.path.basename(fq2).rstrip(".fastq.gz") + "_fastqc", "fastqc_data.txt")


        for folder in [
                output_results_path, 
                qc_path, 
                skewer_path, 
                filtered_reads_path,
                expression_path, 
#                kallisto_path, 
                star_path,
                fusion_path, 
                mapsplice_path, 
#                pizzly_path, 
                fusioncatcher_path, 
                starfusion_path, 
#                starchip_path, 
                infusion_path, 
                soapfuse_path,
                fetchdata_path
            ]:
            IOMethods.create_folder(folder, self.logger)

        # get a list of tools from the samples.db file that have been run previously on this sample
        state_tools = self.samples.get_tool_list_from_state(sample_id)
        # get a list of tools from the config file which shall be run on this sample
        tools = cfg.tools
        cmds = cfg.commands
        module_dir = cfg.module_dir
        # Define cmd strings for each program
        # urla: mapsplice requires gunzip'd read files and process substitutions don't seem to work in slurm scripts...
        #       process substitution do somehow not work from this script - c/p the command line to the terminal, however, works w/o issues?!
        cmd_fastqc = "{0} --nogroup --extract -t 6 -o {1} {2} {3}".format(cmds["fastqc"], qc_path, fq1, fq2)
        cmd_qc_parser = "{0} -i {1} {2} -o {3}".format(os.path.join(module_dir, "misc", "qc_parser.py"), fastqc_1, fastqc_2, qc_table_path)
        cmd_skewer = "{0} -q {1} -i {2} {3} -o {4}".format(os.path.join(module_dir, "tool_wrapper", "skewer_wrapper.py"), qc_table_path, fq1, fq2, skewer_path)

        fq0 = ""
        if "QC" in tools:
            fq0 = os.path.join(skewer_path, "out_file-trimmed.fastq.gz")
            fq1 = os.path.join(skewer_path, "out_file-trimmed-pair1.fastq.gz")
            fq2 = os.path.join(skewer_path, "out_file-trimmed-pair2.fastq.gz")
        else:
            qc_table_path = "None"

        # (0) Readfilter
        cmd_star_filter = "{0} --genomeDir {1} --outFileNamePrefix {2}_ --readFilesCommand zcat --readFilesIn {3} {4} --outFilterMultimapNmax 100 --outSAMmultNmax 1 --chimSegmentMin 10 --chimJunctionOverhangMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax {5} --alignIntronMax {5} --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --seedSearchStartLmax 20 --winAnchorMultimapNmax 50 --outSAMtype BAM Unsorted --chimOutType Junctions WithinBAM --outSAMunmapped Within KeepPairs --runThreadN waiting_for_cpu_number".format(cmds["star"], star_index_path, os.path.join(filtered_reads_path, sample_id), fq1, fq2, cfg.max_dist_proper_pair)
        cmd_read_filter = "{0} --input {1}_Aligned.out.bam --output {1}_Aligned.out.filtered.bam".format(os.path.join(module_dir, "fusionreadfilter.py"), os.path.join(filtered_reads_path, sample_id))
        # re-define fastq's if filtering is on (default)
        fq0 = ""
        if "Readfilter" in tools:
            fq0 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(".fastq.gz", "_filtered_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered.fastq.gz"))

        cmd_bam_to_fastq = "{0} fastq -0 {1} -1 {2} -2 {3} --threads waiting_for_cpu_number {4}_Aligned.out.filtered.bam".format(cmds["samtools"], fq0, fq1, fq2, os.path.join(filtered_reads_path, sample_id))
        # (1) Kallisto expression quantification (required for pizzly)
#        cmd_kallisto = "{0} quant --threads waiting_for_cpu_number --genomebam --gtf {1} --chromosomes {2} --index {3} --fusion --output-dir waiting_for_output_string {4} {5}".format(cmds["kallisto"], genes_gtf_path, genome_sizes_path, kallisto_index_path, fq1, fq2)
        # (2) Star expression quantification (required for starfusion and starchip)
        cmd_star = "{0} --genomeDir {1} --outFileNamePrefix waiting_for_output_string --runThreadN waiting_for_cpu_number --runMode alignReads --readFilesIn {2} {3} --readFilesCommand zcat --chimSegmentMin 10 --chimJunctionOverhangMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax {4} --alignIntronMax {4} --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --seedSearchStartLmax 20 --winAnchorMultimapNmax 50 --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold --chimOutJunctionFormat 1".format(cmds["star"], star_index_path, fq1, fq2, cfg.max_dist_proper_pair)
        # (3) Mapslice
        # urla: the "keep" parameter requires gunzip >= 1.6
        cmd_extr_fastq1 = "gunzip --keep {0}".format(fq1)
        cmd_extr_fastq2 = "gunzip --keep {0}".format(fq2)
        # Added python interpreter to circumvent external hardcoded shell script
        cmd_mapsplice = "{0} --chromosome-dir {1} -x {2} -1 {3} -2 {4} --threads waiting_for_cpu_number --output {5} --qual-scale phred33 --bam --seglen 20 --min-map-len 40 --gene-gtf {6} --fusion".format(cmds["mapsplice"], genome_chrs_path, bowtie_index_path, fq1[:-3], fq2[:-3], mapsplice_path, genes_gtf_path)
        # (4) Fusiocatcher
        cmd_fusioncatcher = "{0} --input {1} --data {2} --output {3} -p waiting_for_cpu_number".format(cmds["fusioncatcher"], ",".join([fq1, fq2]), fusioncatcher_index_path, fusioncatcher_path)
        # star-fusion and star-chip can be run upon a previous star run (this MUST NOT be the star_filter run, but the star_expression run)
        # (5)
        cmd_starfusion = "{0} --chimeric_junction {1} --genome_lib_dir {2} --CPU waiting_for_cpu_number --output_dir {3}".format(cmds["starfusion"], "{}_Chimeric.out.junction".format(os.path.join(star_path, sample_id)), starfusion_index_path, starfusion_path)
        # (7)
#        cmd_starchip = "{0} {1} {2} {3}".format(cmds["starchip"], os.path.join(starchip_path, "starchip"), "{}_Chimeric.out.junction".format(os.path.join(star_path, sample_id)), starchip_param_path)
        # (6) Pizzly
#        cmd_pizzly = "{0} -k 29 --gtf {1} --cache {2} --fasta {3} --output {4} {5}".format(cmds["pizzly"], genes_gtf_path, pizzly_cache_path, genes_fasta_path, os.path.join(pizzly_path, "kallizzy"), os.path.join(kallisto_path, "fusion.txt"))
#        cmd_pizzly2 = "{0} {1} {2}".format(cmds["pizzly_cmd2"], "{}.json".format(os.path.join(pizzly_path, "kallizzy")), "{}.json.txt".format(os.path.join(pizzly_path, "kallizzy")))
        # (8) Infusion
        cmd_infusion = "{0} -1 {1} -2 {2} --skip-finished --min-unique-alignment-rate 0 --min-unique-split-reads 0 --allow-non-coding --out-dir {3} {4}".format(cmds["infusion"], fq1, fq2, infusion_path, infusion_cfg_path)
        # (x) Soapfuse
        cmd_soapfuse = "{0} -q {1} -i {2} -o {3}".format(os.path.join(module_dir, "tool_wrapper", "soapfuse_wrapper.py"), qc_table_path, " ".join([fq1, fq2]), soapfuse_path)
        # (9) Data collection
        cmd_fetchdata = "{0} -i {1} -o {2} -s {3} --fq1 {4} --fq2 {5} --fusion_support {6}".format(os.path.join(module_dir, "fetchdata.py"), output_results_path, fetchdata_path, sample_id, fq1, fq2, tool_num_cutoff)
        # (10) De novo assembly of fusion transcripts
        # urla: This is currently still under active development and has not been tested thoroughly
#        cmd_denovoassembly = "{0} -i waiting_for_gene_list_input -b {1}_Aligned.out.bam -g {2} -t {3} -o waiting_for_assembly_out_dir".format(os.path.join(module_dir, "denovoassembly.py"), os.path.join(filtered_reads_path, sample_id), ref_genome, ref_trans)
        # (X) Sample monitoring
        cmd_samples = "{0} --db_path={1} --sample_id={2} --action=append_state --tool=".format(os.path.join(module_dir, "misc", "samples.py"), self.samples.db_path, sample_id)

        # set final lists of executable tools and path
        exe_tools = [
            "QC", #0
            "Readfilter", #1
#            "Kallisto", #2
            "Star", #3
            "Mapsplice", #4
            "Fusioncatcher", #5
            "Starfusion", #6
#            "Pizzly", #7
#            "Starchip", #8
            "Infusion", #9
            "Soapfuse", #10
            "Fetchdata" #11
#            "Assembly" #12
            ]
        exe_cmds = [
            " && ".join([cmd_fastqc, cmd_qc_parser, cmd_skewer]), #0
            " && ".join([cmd_star_filter, cmd_read_filter, cmd_bam_to_fastq]), #1
#            cmd_kallisto, #2
            cmd_star, #3
            " && ".join([cmd_extr_fastq1, cmd_extr_fastq2, cmd_mapsplice]), #4
            cmd_fusioncatcher, #5
            cmd_starfusion, #6
#            " && ".join([cmd_pizzly, cmd_pizzly2]), #7
#            cmd_starchip, #8
            cmd_infusion, #9
            cmd_soapfuse, #10
            cmd_fetchdata #11
#            cmd_denovoassembly #12
            ]
        exe_path = [
            qc_path, #0
            filtered_reads_path, #1
#            kallisto_path, #2
            star_path, #3
            mapsplice_path, #4
            fusioncatcher_path, #5
            starfusion_path, #6
#            pizzly_path, #7
#            starchip_path, #8
            infusion_path, #9
            soapfuse_path, #10
            fetchdata_path #11
#            "" #12
            ]

        # create and submit slurm job if the tool is requested and hasn't been run before
        for i, tool in enumerate(exe_tools, 0):
            if tool in tools:
                dependency = []
                # check dependencies of the pipeline.
                # Besides tool dependencies (Pizzly -> Kallisto, Starfusion/Starchip -> Star), read filtering is mandatory
                # Processing will be skipped if a certain dependency was not found (either pre-processed data of the configs tool string are checked)
                if tool in state_tools:
                    # urla: the primary idea behind this flag is to allow multiple fetchdata executions during processing
                    #       nevertheless, re-processing of the same data with a newer version of a tool will also be straightforward (but overwriting previous results, of course)
#                    if self.overwrite:
#                        self.logger.info("Executing {0} although it looks like a previous run finished successfully. Results in {1} may be overwritten".format(tool, exe_path[i]))
#                    else:
                    self.logger.info("Skipping {0} as it looks like a previous run finished successfully. Results should be in {1}".format(tool, exe_path[i]))
                    continue
                else:
                    if tool == "Readfilter" and "Readfilter" not in tools:
                        self.logger.error(
                                """Error 99: Sample {} will be skipped due to missing read filtering.\n
                                Read filtering is currently a mandatory step for the processing.\n
                                Because you haven't run it before for this sample, you have to include \"Readfilter\" in the tool selection in your config.\n
                                """.format(sample_id))
                        print("Error 99: Sample {} will be skipped due to missing read filtering.".format(sample_id))
                        return 0
                    elif tool == "Pizzly" and "Kallisto" not in tools:
                        self.logger.error(
                                """Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.\n
                                Pizzly builds on Kallisto and it is therefore mandatory to run this first.\n
                                Because you haven't run it before for this sample, you have to include \"Kallisto\" in the tool selection in your config.\n
                                """.format(sample_id))
                        print("Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.".format(tool, sample_id))
                        continue
                    elif (tool == "Starfusion" or tool == "Starchip") and "Star" not in tools:
                        self.logger.error(
                                """Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.\n
                                {0} builds on Star and it is therefore mandatory to run this first.\n
                                Because you haven't run it before for this sample, you have to include \"Star\" in the tool selection in your config.\n
                                """.format(sample_id))
                        print("Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.".format(tool, sample_id))
                        continue

                # prepare slurm jobs: get ressources, create uid, set output path and check dependencies
                self.logger.debug("Submitting {} run to slurm".format(tool))
                cpu = cfg.resources[tool.lower()]["cpu"]
                mem = cfg.resources[tool.lower()]["mem"]
                uid = "-".join([cfg.pipeline_name, tool, sample_id])
                if tool == "Star":
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string", "{}_".format(os.path.join(exe_path[i], sample_id))).replace("waiting_for_cpu_number", str(cpu))
                else:
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string", exe_path[i]).replace("waiting_for_cpu_number", str(cpu))
                cmd = " && ".join([exe_cmds[i], cmd_samples + tool])
                # Managing slurm dependencies
                que_sys = cfg.queueing_system
                if tool == "Pizzly":
                    dependency = Queueing.get_jobs_by_name("Kallisto-{0}".format(sample_id), que_sys)
                elif tool == "Starfusion" or tool == "Starchip":
                    dependency = Queueing.get_jobs_by_name("Star-{0}".format(sample_id), que_sys)
                elif tool == "Fetchdata":
                    dependency = Queueing.get_jobs_by_name(sample_id, que_sys)
                elif tool == "Assembly":
                    dependency = Queueing.get_jobs_by_name("Fetchdata-{0}".format(sample_id), que_sys)
                elif tool == "ReadFilter":
                    dependency = Queueing.get_jobs_by_name("QC-{0}".format(sample_id), que_sys)
                dependency.extend(Queueing.get_jobs_by_name("Readfilter-{0}".format(sample_id), que_sys))
                dependency.extend(Queueing.get_jobs_by_name("QC-{0}".format(sample_id), que_sys))
                self.logger.debug("Submitting slurm job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd, exe_path[i], dependency))
                self.submit_job(uid, cmd, cpu, mem, exe_path[i], dependency, "")
            else:
                self.logger.info("Skipping {0} as it is not selected for execution (Selected are: {1})".format(tool, tools))

    def submit_job(self, uid, cmd, cores, mem_usage, output_results_folder, dependencies, mail):
        """Submit job to slurm scheduling"""
        que_sys = cfg.queueing_system
        already_running = Queueing.get_jobs_by_name(uid, que_sys)
        if not already_running:
            # urla: for compatibility reasons (and to be independent of shell commands), concatenated commands are splitted again,
            #       dependencies within the splitted groups updated and everything submitted sequentially to the queueing system
            module_file = os.path.join(cfg.module_dir, "build_env.sh")

            for i, cmd_split in enumerate(cmd.split(" && ")):
                if not que_sys in ["slurm", "pbs"]:
                    cmd_split = cmd_split.split(" ")
                dependencies.extend(Queueing.get_jobs_by_name("{0}_CMD{1}".format(uid, i - 1), que_sys))
                Queueing.submit("{0}_CMD{1}".format(uid, i), cmd_split, cores, mem_usage, output_results_folder, dependencies, cfg.partition, cfg.user, cfg.time_limit, mail, module_file, que_sys)
                time.sleep(0.5)
        else:
            self.logger.error("A job with this application/sample combination is currently running. Skipping {} in order to avoid unintended data loss.".format(uid))

def main():
    """Parse command line arguments and start script"""
    parser = argparse.ArgumentParser(description='Processing of demultiplexed FASTQs')
    # required arguments
    parser.add_argument('-i', '--input', dest='input_paths', nargs='+', help='Specify full path of the fastq folder to process.', required=True)
    parser.add_argument('-o', '--output-folder', dest='output_folder', help='Specify full path of the folder to save the results into.', required=True)
    # optional arguments
    parser.add_argument('--tool_support', dest='tool_support', help='The number of tools which are required to support a fusion event.', default="1")
    parser.add_argument('--version', dest='version', help='Get version info')
    args = parser.parse_args()

    # if version is request, print it and exit
    if args.version:
        print(cfg.__version__)
        sys.exit(0)

    script_call = "python {} -i {} -o {}".format(os.path.realpath(__file__), " ".join([os.path.abspath(x) for x in args.input_paths]), os.path.abspath(args.output_folder))

    proc = Processing(script_call, args.input_paths, args.output_folder)
    proc.run(args.tool_support)

    with open(os.path.join(args.output_folder, "process.sh"), "w") as outf:
        outf.write("#!/bin/sh\n\n{}".format(script_call))

if __name__ == '__main__':
    main()
