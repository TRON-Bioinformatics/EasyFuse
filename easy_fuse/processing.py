#!/usr/bin/env python3

"""
Read fastq files
create per-sample (sub)folder tree
run gene fusion predictions
start data collection
during all steps, perform slurm scheduling and process monitoring

@author: Tron (PASO), BNT (URLA)
"""
import os
import os.path
import time
from argparse import ArgumentParser
from shutil import copy

import logzero
from logzero import logger

import easy_fuse.misc.io_methods as io_methods
from easy_fuse import __version__
from easy_fuse.misc import queueing
from easy_fuse.misc.config import EasyFuseConfiguration
from easy_fuse.misc.samples import SamplesDB


class Processing(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""

    def __init__(self, cmd, input_paths, working_dir, config: EasyFuseConfiguration, jobname_suffix):
        """Parameter initiation and work folder creation. Start of progress logging."""
        self.working_dir = os.path.abspath(working_dir)
        logzero.logfile(os.path.join(self.working_dir, "easyfuse_processing.log"))
        io_methods.create_folder(self.working_dir)

        logger.info("Starting easyfuse: CMD - {}".format(cmd))
        self.input_paths = [os.path.abspath(file) for file in input_paths]
        self.samples = SamplesDB(os.path.join(self.working_dir, "samples.db"))

        self.cfg = config
        copy(self.cfg.config_file, working_dir)

        self.jobname_suffix = None
        if jobname_suffix:
            self.jobname_suffix = jobname_suffix


    # The run method simply greps and organises fastq input files.
    # Fastq pairs (single end input is currently not supported) are then send to "execute_pipeline"
    def run(self, tool_num_cutoff):
        """General parameter setting, identification of fastq files and initiation of processing"""
        logger.info("Pipeline Version: {}".format(__version__))

        ref_genome = self.cfg["general"]["ref_genome_build"]
        ref_trans = self.cfg["general"]["ref_trans_version"]

        logger.info("Reference Genome: {0}, Reference Transcriptome: {1}".format(ref_genome, ref_trans))

        # get fastq files
        fastqs = io_methods.get_fastq_files(self.input_paths)
        sample_list = io_methods.pair_fastq_files(fastqs)
        for sample in sample_list:
            if sample[1] and sample[2]:
                logger.info("Processing Sample ID: {} (paired end)".format(sample[0]))
                logger.info("Sample 1: {}".format(sample[1]))
                logger.info("Sample 2: {}".format(sample[2]))
                self.execute_pipeline(sample[1], sample[2], sample[0], ref_genome, ref_trans, tool_num_cutoff)

        # summarize all data if selected
        if "Summary" in self.cfg["general"]["tools"].split(","):
            dependency = []
            for sample in sample_list:
                dependency.extend(
                    queueing.get_jobs_by_name("Fetchdata-{}".format(sample[0]), self.cfg["general"]["queueing_system"]))
            modelling_string = ""
            if self.cfg["other_files"]["easyfuse_model"]:
                modelling_string = " --model_predictions"

            cmd_summarize = "summarize_data --input {0}{1} -c {2}".format(
                self.working_dir,
                modelling_string,
                self.cfg.config_file)

            logger.debug(
                "Submitting job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd_summarize, self.working_dir, dependency))
            resources_summary = self.cfg["resources"]["summary"].split(",")
            cpu = resources_summary[0]
            mem = resources_summary[1]
            self.submit_job("-".join([self.cfg["general"]["pipeline_name"], "Summary", str(int(round(time.time())))]),
                            cmd_summarize, cpu, mem, self.working_dir, dependency, self.cfg["general"]["receiver"])

    # Per sample, define input parameters and execution commands, create a folder tree and submit runs to slurm
    def execute_pipeline(self, fq1, fq2, sample_id, ref_genome, ref_trans, tool_num_cutoff):
        """Create sample specific subfolder structure and run tools on fastq files"""
        self.samples.add_sample(sample_id, "NA", fq1, fq2)

        refs = self.cfg["references"]

        # Genome/Gene references to use
        genome_chrs_path = refs["genome_fastadir"]
        genes_gtf_path = refs["genes_gtf"]

        # Path' to specific indices
        indices = self.cfg["indices"]
        other_files = self.cfg["other_files"]

        bowtie_index_path = indices["bowtie"]
        star_index_path = indices["star"]
        starfusion_index_path = indices["starfusion"]
        fusioncatcher_index_path = indices["fusioncatcher"]
        infusion_cfg_path = other_files["infusion_cfg"]

        # Output results folder creation - currently included:
        # 1) Gene/Isoform expression: star
        # 2) Fusion prediction: mapsplice, fusioncatcher, starfusion, infusion, soapfuse
        output_results_path = os.path.join(self.working_dir, "Sample_{}".format(sample_id))
        qc_path = os.path.join(output_results_path, "qc")
        skewer_path = os.path.join(qc_path, "skewer")
        qc_table_path = os.path.join(qc_path, "qc_table.txt")
        filtered_reads_path = os.path.join(output_results_path, "filtered_reads")
        expression_path = os.path.join(output_results_path, "expression")
        star_path = os.path.join(expression_path, "star")
        fusion_path = os.path.join(output_results_path, "fusion")
        mapsplice_path = os.path.join(fusion_path, "mapsplice")
        fusioncatcher_path = os.path.join(fusion_path, "fusioncatcher")
        starfusion_path = os.path.join(fusion_path, "starfusion")
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
            star_path,
            fusion_path,
            mapsplice_path,
            fusioncatcher_path,
            starfusion_path,
            infusion_path,
            soapfuse_path,
            fetchdata_path
        ]:
            io_methods.create_folder(folder)

        # get a list of tools from the samples.db file that have been run previously on this sample
        state_tools = self.samples.get_tool_list_from_state(sample_id)
        # get a list of tools from the config file which shall be run on this sample
        tools = self.cfg["general"]["tools"].split(",")
        cmds = self.cfg["commands"]
        # Define cmd strings for each program
        # urla: mapsplice requires gunzip'd read files and process substitutions don't seem to work in slurm scripts...
        #       process substitution do somehow not work from this script - c/p the command line to the terminal, however, works w/o issues?!
        # (0) QC
        cmd_fastqc = "{0} --nogroup --extract -t 6 -o {1} {2} {3}".format(cmds["fastqc"], qc_path, fq1, fq2)
        cmd_qc_parser = "qc_parser " \
                        "-i {0} {1} " \
                        "-o {2}".format(
            fastqc_1,
            fastqc_2,
            qc_table_path)

        cmd_skewer = "skewer_wrapper -q {0} -i {1} {2} -o {3} -b {4} -m {5}".format(
            qc_table_path, fq1, fq2, skewer_path,
            self.cfg["commands"]["skewer"], self.cfg["general"]["min_read_len_perc"])

        fq0 = ""
        if "QC" in tools:
            fq0 = os.path.join(skewer_path, "out_file-trimmed.fastq.gz")
            fq1 = os.path.join(skewer_path, "out_file-trimmed-pair1.fastq.gz")
            fq2 = os.path.join(skewer_path, "out_file-trimmed-pair2.fastq.gz")
        else:
            qc_table_path = "None"

        # (1) Readfilter
        cmd_star_filter = "{0} " \
                          "--genomeDir {1} " \
                          "--outFileNamePrefix {2}_ " \
                          "--readFilesCommand zcat " \
                          "--readFilesIn {3} {4} " \
                          "--outFilterMultimapNmax 100 " \
                          "--outSAMmultNmax 1 " \
                          "--chimSegmentMin 10 " \
                          "--chimJunctionOverhangMin 10 " \
                          "--alignSJDBoverhangMin 10 " \
                          "--alignMatesGapMax {5} " \
                          "--alignIntronMax {5} " \
                          "--chimSegmentReadGapMax 3 " \
                          "--alignSJstitchMismatchNmax 5 -1 5 5 " \
                          "--seedSearchStartLmax 20 " \
                          "--winAnchorMultimapNmax 50 " \
                          "--outSAMtype BAM Unsorted " \
                          "--chimOutType Junctions WithinBAM " \
                          "--outSAMunmapped Within KeepPairs " \
                          "--runThreadN waiting_for_cpu_number".format(
            cmds["star"], star_index_path, os.path.join(filtered_reads_path, sample_id), fq1, fq2,
            self.cfg["general"]["max_dist_proper_pair"])

        cmd_read_filter = "fusionreadfilter --input {0}_Aligned.out.bam --output {0}_Aligned.out.filtered.bam".format(
            os.path.join(filtered_reads_path, sample_id))

        # re-define fastq's if filtering is on (default)
        fq0 = ""
        if "Readfilter" in tools:
            fq0 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(".fastq.gz",
                                                                                                      "_filtered_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered.fastq.gz"))

        cmd_bam_to_fastq = "{0} fastq -0 {1} -1 {2} -2 {3} --threads waiting_for_cpu_number {4}_Aligned.out.filtered.bam".format(
            cmds["samtools"], fq0, fq1, fq2, os.path.join(filtered_reads_path, sample_id))

        # (2) Star expression quantification (required for starfusion)
        cmd_star = "{0} --genomeDir {1} --outFileNamePrefix waiting_for_output_string --runThreadN waiting_for_cpu_number --runMode alignReads --readFilesIn {2} {3} --readFilesCommand zcat --chimSegmentMin 10 --chimJunctionOverhangMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax {4} --alignIntronMax {4} --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --seedSearchStartLmax 20 --winAnchorMultimapNmax 50 --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold --chimOutJunctionFormat 1".format(
            cmds["star"], star_index_path, fq1, fq2, self.cfg["general"]["max_dist_proper_pair"])

        # (3) Mapslice
        # urla: the "keep" parameter requires gunzip >= 1.6
        cmd_extr_fastq1 = "gunzip --keep {0}".format(fq1)
        cmd_extr_fastq2 = "gunzip --keep {0}".format(fq2)

        # Added python interpreter to circumvent external hardcoded shell script
        cmd_mapsplice = "python2 {0} --chromosome-dir {1} -x {2} -1 {3} -2 {4} --threads waiting_for_cpu_number --output {5} --qual-scale phred33 --bam --seglen 20 --min-map-len 40 --gene-gtf {6} --fusion".format(
            cmds["mapsplice"], genome_chrs_path, bowtie_index_path, fq1[:-3], fq2[:-3], mapsplice_path, genes_gtf_path)

        # (4) Fusioncatcher
        cmd_fusioncatcher = "{0} --input {1} --data {2} --output {3} -p waiting_for_cpu_number".format(
            cmds["fusioncatcher"], ",".join([fq1, fq2]), fusioncatcher_index_path, fusioncatcher_path)

        # (5) Starfusion
        cmd_starfusion = "{0} --chimeric_junction {1} --genome_lib_dir {2} --CPU waiting_for_cpu_number --output_dir {3}".format(
            cmds["starfusion"], "{}_Chimeric.out.junction".format(os.path.join(star_path, sample_id)),
            starfusion_index_path, starfusion_path)

        # (6) Infusion
        cmd_infusion = "{0} -1 {1} -2 {2} --skip-finished --min-unique-alignment-rate 0 --min-unique-split-reads 0 --allow-non-coding --out-dir {3} {4}".format(
            cmds["infusion"], fq1, fq2, infusion_path, infusion_cfg_path)

        # (7) Soapfuse
        cmd_soapfuse = "soapfuse_wrapper -q {0} -i {1} -o {2} -b {3} -c {4}".format(
            qc_table_path, " ".join([fq1, fq2]),
            soapfuse_path, self.cfg["commands"]["soapfuse"], self.cfg["other_files"]["soapfuse_cfg"])

        # (8) Data collection
        cmd_fetchdata = "fetchdata -i {0} -o {1} -s {2} -c {3} --fq1 {4} --fq2 {5} --fusion_support {6}".format(
            output_results_path, fetchdata_path, sample_id, self.cfg.config_file, fq1, fq2, tool_num_cutoff)

        # (X) Sample monitoring
        cmd_samples = "samples --db_path={0} --sample_id={1} --action=append_state --tool=".format(
            self.samples.db_path, sample_id)

        # set final lists of executable tools and path
        exe_tools = [
            "QC",  # 0
            "Readfilter",  # 1
            "Star",  # 2
            "Mapsplice",  # 3
            "Fusioncatcher",  # 4
            "Starfusion",  # 5
            "Infusion",  # 6
            "Soapfuse",  # 7
            "Fetchdata"  # 8
        ]
        exe_cmds = [
            " && ".join([cmd_fastqc, cmd_qc_parser, cmd_skewer]),  # 0
            " && ".join([cmd_star_filter, cmd_read_filter, cmd_bam_to_fastq]),  # 1
            cmd_star,  # 2
            " && ".join([cmd_extr_fastq1, cmd_extr_fastq2, cmd_mapsplice]),  # 3
            cmd_fusioncatcher,  # 4
            cmd_starfusion,  # 5
            cmd_infusion,  # 6
            cmd_soapfuse,  # 7
            cmd_fetchdata  # 8
        ]
        exe_path = [
            qc_path,  # 0
            filtered_reads_path,  # 1
            star_path,  # 2
            mapsplice_path,  # 3
            fusioncatcher_path,  # 4
            starfusion_path,  # 5
            infusion_path,  # 6
            soapfuse_path,  # 7
            fetchdata_path  # 8
        ]

        # create and submit slurm job if the tool is requested and hasn't been run before
        for i, tool in enumerate(exe_tools, 0):
            if tool in tools:
                dependency = []
                # check dependencies of the pipeline.
                # Besides tool dependencies (Pizzly -> Kallisto, Starfusion/Starchip -> Star), read filtering is mandatory
                # Processing will be skipped if a certain dependency was not found (either pre-processed data of the configs tool string are checked)
                if tool in state_tools:
                    logger.info(
                        "Skipping {0} as it looks like a previous run finished successfully. Results should be in {1}".format(
                            tool, exe_path[i]))
                    continue
                else:
                    if tool == "Readfilter" and "Readfilter" not in tools:
                        logger.error(
                            """Error 99: Sample {} will be skipped due to missing read filtering.\n
                            Read filtering is currently a mandatory step for the processing.\n
                            Because you haven't run it before for this sample, you have to include \"Readfilter\" in the tool selection in your config.\n
                            """.format(sample_id))
                        print("Error 99: Sample {} will be skipped due to missing read filtering.".format(sample_id))
                        return 0
                    elif tool == "Pizzly" and "Kallisto" not in tools:
                        logger.error(
                            """Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.\n
                            Pizzly builds on Kallisto and it is therefore mandatory to run this first.\n
                            Because you haven't run it before for this sample, you have to include \"Kallisto\" in the tool selection in your config.\n
                            """.format(sample_id))
                        print(
                            "Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.".format(
                                tool, sample_id))
                        continue
                    elif (tool == "Starfusion" or tool == "Starchip") and "Star" not in tools:
                        logger.error(
                            """Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.\n
                            {0} builds on Star and it is therefore mandatory to run this first.\n
                            Because you haven't run it before for this sample, you have to include \"Star\" in the tool selection in your config.\n
                            """.format(sample_id))
                        print(
                            "Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.".format(
                                tool, sample_id))
                        continue

                # prepare slurm jobs: get ressources, create uid, set output path and check dependencies
                logger.debug("Submitting {} run".format(tool))
                cpumem = self.cfg["resources"][tool.lower()].split(",")
                cpu = cpumem[0]
                mem = cpumem[1]
                uid = "-".join([self.cfg["general"]["pipeline_name"], tool, sample_id])
                if tool == "Star":
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string",
                                                      "{}_".format(os.path.join(exe_path[i], sample_id))).replace(
                        "waiting_for_cpu_number", str(cpu))
                else:
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string", exe_path[i]).replace(
                        "waiting_for_cpu_number", str(cpu))
                cmd = " && ".join([exe_cmds[i], cmd_samples + tool])
                # Managing slurm dependencies
                que_sys = self.cfg["general"]["queueing_system"]
                if tool == "Pizzly":
                    dependency = queueing.get_jobs_by_name("Kallisto-{0}".format(sample_id), que_sys)
                elif tool == "Starfusion" or tool == "Starchip":
                    dependency = queueing.get_jobs_by_name("Star-{0}".format(sample_id), que_sys)
                elif tool == "Fetchdata":
                    dependency = queueing.get_jobs_by_name(sample_id, que_sys)
                elif tool == "Assembly":
                    dependency = queueing.get_jobs_by_name("Fetchdata-{0}".format(sample_id), que_sys)
                elif tool == "ReadFilter":
                    dependency = queueing.get_jobs_by_name("QC-{0}".format(sample_id), que_sys)
                dependency.extend(queueing.get_jobs_by_name("Readfilter-{0}".format(sample_id), que_sys))
                dependency.extend(queueing.get_jobs_by_name("QC-{0}".format(sample_id), que_sys))
                logger.debug("Submitting job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd, exe_path[i], dependency))
                self.submit_job(uid, cmd, cpu, mem, exe_path[i], dependency, "")
            else:
                logger.info("Skipping {0} as it is not selected for execution (Selected are: {1})".format(tool, tools))

    def submit_job(self, uid, cmd, cores, mem_usage, output_results_folder, dependencies, mail):
        """Submit job to for process generation"""
        que_sys = self.cfg["general"]["queueing_system"]
        already_running = queueing.get_jobs_by_name(uid, que_sys)
        if not already_running:
            # urla: for compatibility reasons (and to be independent of shell commands), concatenated commands are splitted again,
            #       dependencies within the splitted groups updated and everything submitted sequentially to the queueing system
            module_file = self.cfg["general"]["build_env"]

            for i, cmd_split in enumerate(cmd.split(" && ")):
                if not que_sys in ["slurm", "pbs"]:
                    cmd_split = cmd_split.split(" ")
                dependencies.extend(queueing.get_jobs_by_name("{0}_CMD{1}".format(uid, i - 1), que_sys))
                if self.jobname_suffix:
                    queueing.submit("{0}_CMD{1}-{2}".format(uid, i, self.jobname_suffix), cmd_split, cores, mem_usage,
                                    output_results_folder, dependencies, self.cfg["general"]["partition"],
                                    self.cfg["general"]["user"], self.cfg["general"]["time_limit"], mail, module_file,
                                    que_sys)
                else:
                    queueing.submit("{0}_CMD{1}".format(uid, i), cmd_split, cores, mem_usage, output_results_folder,
                                    dependencies, self.cfg["general"]["partition"], self.cfg["general"]["user"],
                                    self.cfg["general"]["time_limit"], mail, module_file, que_sys)
                time.sleep(0.5)
        else:
            logger.error(
                "A job with this application/sample combination is currently running. Skipping {} in order to avoid unintended data loss.".format(
                    uid))


def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(prog='EasyFuse', description='Run EasyFuse pipeline')
    # required arguments
    parser.add_argument('-i', '--input', dest='input_paths', nargs='+',
                        help='Specify full path of the fastq folder to process.', required=True)
    parser.add_argument('-o', '--output-folder', dest='output_folder',
                        help='Specify full path of the folder to save the results into.', required=True)
    # optional arguments
    parser.add_argument('-c', '--config-file', dest='config_file', required=True,
                        help='Specify alternative config file to use for your analysis')
    parser.add_argument('--tool_support', dest='tool_support',
                        help='The number of tools which are required to support a fusion event.', default="1")
    parser.add_argument('-p', '--jobname_suffix', dest='jobname_suffix',
                        help='Specify a jobname suffix for running jobs on a queueing system', default='')
    parser.add_argument('-V', '--version', dest='version', action='version',
                        version="%(prog)s {version}".format(version=__version__), help='Get version info')
    args = parser.parse_args()

    jobname_suffix = ""
    if args.jobname_suffix:
        jobname_suffix = "-p {}".format(args.jobname_suffix)

    config = EasyFuseConfiguration(args.config_file)

    # records CLI call
    # TODO: think of better ways of recording what the content of files (ie: config, FASTQs, etc.) was as files do
    #  change and only paths are recorded here
    script_call = "easy-fuse -i {} {} {} -o {}".format(
        " ".join([os.path.abspath(x) for x in args.input_paths]),
        "-c {}".format(config.config_file),
        jobname_suffix,
        os.path.abspath(args.output_folder))
    with open(os.path.join(args.output_folder, "process.sh"), "w") as outf:
        outf.write("#!/bin/sh\n\n{}".format(script_call))

    proc = Processing(script_call, args.input_paths, args.output_folder, config, args.jobname_suffix)
    proc.run(args.tool_support)



# if __name__ == '__main__':
#    main()
