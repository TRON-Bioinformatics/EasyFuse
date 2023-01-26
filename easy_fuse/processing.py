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
from shutil import copy

import logzero
from logzero import logger

import easy_fuse
import easy_fuse.misc.io_methods as io_methods
from easy_fuse import __version__
from easy_fuse.misc import queueing
from easy_fuse.misc.config import EasyFuseConfiguration

DEFAULT_CPU_MEM = "1,16"

READ_FILTER_STEP = "readfilter"


class Processing(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""

    def __init__(self, cmd, input_paths, working_dir, config: EasyFuseConfiguration, jobname_suffix):
        """Parameter initiation and work folder creation. Start of progress logging."""
        self.working_dir = os.path.abspath(working_dir)
        self.log_path = os.path.join(self.working_dir, "easyfuse_processing.log")
        logzero.logfile(self.log_path)
        io_methods.create_folder(self.working_dir)

        logger.info("Starting easyfuse: CMD - {}".format(cmd))
        self.input_paths = [os.path.abspath(file) for file in input_paths]

        self.cfg = config
        copy(self.cfg.config_file, working_dir)

        self.jobname_suffix = None
        if jobname_suffix:
            self.jobname_suffix = jobname_suffix


    # The run method simply greps and organises fastq input files.
    # Fastq pairs (single end input is currently not supported) are then send to "execute_pipeline"
    def run(self):
        """General parameter setting, identification of fastq files and initiation of processing"""
        logger.info("Pipeline Version: {}".format(__version__))

        ref_genome = self.cfg["general"]["ref_genome_build"]
        ref_trans = self.cfg["general"]["ref_trans_version"]

        logger.info("Reference Genome: {0}, Reference Transcriptome: {1}".format(ref_genome, ref_trans))

        # get fastq files
        fastqs = io_methods.get_fastq_files(self.input_paths)
        sample_list = io_methods.pair_fastq_files(fastqs)
        for sample_id, fq1, fq2 in sample_list:
            if fq1 and fq2:
                logger.info("Processing Sample ID: {} (paired end)".format(sample_id))
                logger.info("Sample 1: {}".format(fq1))
                logger.info("Sample 2: {}".format(fq2))
                self.execute_pipeline(fq1, fq2, sample_id, ref_genome, ref_trans)

        # summarize all data if selected
        if "summary" in self.cfg["general"]["tools"].split(","):
            dependency = []
            for sample in sample_list:
                dependency.extend(
                    queueing.get_jobs_by_name("requantifyBest-{}".format(sample[0]), self.cfg["general"]["queueing_system"]))

            # TODO: add support for running without model and multiple models
            modelling_string = " --model_predictions"

            cmd_summarize = "easy-fuse summarize-data --input {0}{1} -c {2} --samples {3}".format(
                self.working_dir,
                modelling_string,
                self.cfg.config_file,
                " ".join([sample_id for sample_id, _, _ in sample_list])
            )

            logger.debug(
                "Submitting job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd_summarize, self.working_dir, dependency))
            resources_summary = self.cfg["resources"]["summary"].split(",")
            cpu = resources_summary[0]
            mem = resources_summary[1]
            self.submit_job("-".join([self.cfg["general"]["pipeline_name"], "Summary", str(int(round(time.time())))]),
                            cmd_summarize, cpu, mem, self.working_dir, dependency, self.cfg["general"]["receiver"])

    # Per sample, define input parameters and execution commands, create a folder tree and submit runs to slurm
    def execute_pipeline(self, fq1, fq2, sample_id, ref_genome, ref_trans):
        """Create sample specific subfolder structure and run tools on fastq files"""

        refs = self.cfg["references"]

        # Genome/Gene references to use
        genome_fasta = refs["genome_fasta"]
        genome_chrs_path = refs["genome_fastadir"]
        genes_gtf_path = refs["genes_gtf"]
        genes_tsl_path = refs["genes_tsl"]
        genes_adb_path = refs["genes_adb"]


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
        fetchdata_path = os.path.join(output_results_path, "fetchdata")
        # former fetchdata.py
        detected_fusions_path = os.path.join(fetchdata_path, "fetched_fusions")
        detected_fusions_file = os.path.join(detected_fusions_path, "Detected_Fusions.csv")
        context_seq_path = os.path.join(fetchdata_path, "fetched_contextseqs")
        context_seq_file = os.path.join(context_seq_path, "Context_Seqs.csv")
        context_seq_fasta = os.path.join(context_seq_path, "Context_Seqs.csv.fasta")
        star_genome_path = os.path.join(context_seq_path, "STAR_idx")
        star_align_path = os.path.join(context_seq_path, "STAR_align")
        star_align_file = os.path.join(star_align_path, "{}_".format(sample_id))
        classification_path = os.path.join(fetchdata_path, "classification")
        classification_file = os.path.join(classification_path, "classification")

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
            fetchdata_path,
            detected_fusions_path,
            context_seq_path,
            star_genome_path,
            star_align_path,
            classification_path
        ]:
            io_methods.create_folder(folder)

        # get a list of tools from the config file which shall be run on this sample
        tools = self.cfg["general"]["tools"].split(",")
        cmds = self.cfg["commands"]
        # Define cmd strings for each program
        # urla: mapsplice requires gunzip'd read files and process substitutions don't seem to work in slurm scripts...
        #       process substitution do somehow not work from this script - c/p the command line to the terminal, however, works w/o issues?!
        # (0) QC
        cmd_fastqc = "{0} --nogroup --extract -t 6 -o {1} {2} {3}".format(cmds["fastqc"], qc_path, fq1, fq2)
        cmd_qc_parser = "easy-fuse qc-parser " \
                        "-i {0} {1} " \
                        "-o {2}".format(
            fastqc_1,
            fastqc_2,
            qc_table_path)

        cmd_skewer = "easy-fuse skewer-wrapper -q {0} -i {1} {2} -o {3} -b {4} -m {5}".format(
            qc_table_path, fq1, fq2, skewer_path,
            self.cfg["commands"]["skewer"], self.cfg["general"]["min_read_len_perc"])

        fq0 = ""    # TODO: change by None or get rid of
        if "qc" in tools:
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

        cmd_read_filter = "easy-fuse read-filter --input {0}_Aligned.out.bam --output {0}_Aligned.out.filtered.bam".format(
            os.path.join(filtered_reads_path, sample_id))

        # re-define fastq's if filtering is on (default)
        fq0 = ""
        if READ_FILTER_STEP in tools:
            fq0 = os.path.join(
                filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(
                    ".fastq.gz", "_filtered_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered.fastq.gz"))

        cmd_bam_to_fastq = "{0} fastq -0 {1} -1 {2} -2 {3} --threads waiting_for_cpu_number {4}_Aligned.out.filtered.bam".format(
            cmds["samtools"], fq0, fq1, fq2, os.path.join(filtered_reads_path, sample_id))

        # (2) Star expression quantification (required for starfusion)
        cmd_star = "{0} --genomeDir {1} --outFileNamePrefix waiting_for_output_string --runThreadN waiting_for_cpu_number --runMode alignReads --readFilesIn {2} {3} --readFilesCommand zcat --chimSegmentMin 10 --chimJunctionOverhangMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax {4} --alignIntronMax {4} --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --seedSearchStartLmax 20 --winAnchorMultimapNmax 50 --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold --chimOutJunctionFormat 1".format(
            cmds["star"], star_index_path, fq1, fq2, self.cfg["general"]["max_dist_proper_pair"])

        # (3) Mapslice
        # Using "gunzip -c" instead of "gunzip --keep" to keep compatibility with gunzip versions
        cmd_extr_fastq1 = "gunzip -c -f {0} > {1}".format(fq1, fq1[:-3])
        cmd_extr_fastq2 = "gunzip -c -f {0} > {1}".format(fq2, fq2[:-3])

        # Added python interpreter to circumvent external hardcoded shell script
        cmd_mapsplice = "{0} --chromosome-dir {1} -x {2} -1 {3} -2 {4} --threads waiting_for_cpu_number --output {5} --qual-scale phred33 --bam --seglen 20 --min-map-len 40 --gene-gtf {6} --fusion".format(cmds["mapsplice"], genome_chrs_path, bowtie_index_path, fq1[:-3], fq2[:-3], mapsplice_path, genes_gtf_path)
        # (4) Fusioncatcher
        cmd_fusioncatcher = "{0} --input {1} --data {2} --output {3} -p waiting_for_cpu_number".format(
            cmds["fusioncatcher"], ",".join([fq1, fq2]), fusioncatcher_index_path, fusioncatcher_path)

        # (5) Starfusion
        cmd_starfusion = "{0} --chimeric_junction {1} --genome_lib_dir {2} --CPU waiting_for_cpu_number --output_dir {3}".format(
            cmds["starfusion"], "{}_Chimeric.out.junction".format(os.path.join(star_path, sample_id)),
            starfusion_index_path, starfusion_path)

        # (6) Infusion
        # TODO: Use these command line args for more sensitive results:
        # --min-split-reads 0 --min-fragments 0 --min-not-rescued 2
        cmd_infusion = "{0} {1} -1 {2} -2 {3} --skip-finished --min-unique-alignment-rate 0 " \
                       "--min-unique-split-reads 0 --allow-non-coding --out-dir {4} {5}".format(
            cmds["python2"],
            cmds["infusion"],
            fq1,
            fq2,
            infusion_path,
            infusion_cfg_path)

        # (7) Soapfuse
        cmd_soapfuse = "easy-fuse soapfuse-wrapper -q {0} -i {1} -o {2} -b {3} -c {4}".format(
            qc_table_path, " ".join([fq1, fq2]),
            soapfuse_path, self.cfg["commands"]["soapfuse"],
            self.cfg["other_files"]["soapfuse_cfg"])

        # (8) Data collection
        # TODO: embed this as operations within the package
        cmd_readcounts = "easy-fuse count-reads " \
                         "-i {0} " \
                         "-f {1} " \
                         "-o {2}".format(
            os.path.join(filtered_reads_path, "{}_Log.final.out".format(sample_id)),
            "star",
            os.path.join(classification_path, "Star_org_input_reads.txt"))

        cmd_fusiondata = "easy-fuse fusion-parser " \
                         "-i {0} " \
                         "-o {1} " \
                         "-s {2} " \
                         "-f {3} " \
                         "-l {4}".format(
            output_results_path,
            detected_fusions_path,
            sample_id,
            self.cfg["general"]["fusiontools"],
            self.log_path)

        cmd_contextseq = "easy-fuse annotation " \
                         "--detected_fusions {0} " \
                         "--annotation_db {1} " \
                         "--out_csv {2} " \
                         "--genome_fasta {3} " \
                         "--tsl_info {4} " \
                         "--cis_near_dist {5} " \
                         "--context_seq_len {6} " \
                         "--tsl_filter_level {7}".format(
            detected_fusions_file,
            genes_adb_path,
            context_seq_file,
            genome_fasta,
            genes_tsl_path,
            self.cfg["general"]["cis_near_distance"],
            self.cfg["general"]["context_seq_len"],
            self.cfg["general"]["tsl_filter"])

        cmd_starindex = "easy-fuse star-index -i {0} -o {1} -t waiting_for_cpu_number -b {2}".format(
            context_seq_fasta,
            star_genome_path,
            cmds["star"])

        cmd_staralign_fltr = "{0} --genomeDir {1} --readFilesCommand zcat --readFilesIn {2} {3} " \
                             "--outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax -1 " \
                             "--outSAMattributes Standard --outSAMunmapped None " \
                             "--outFilterMismatchNoverLmax 0.02 --runThreadN waiting_for_cpu_number " \
                             "--outFileNamePrefix {4}fltr_ --limitBAMsortRAM 48000000000".format(
            cmds["star"],
            star_genome_path,
            fq1,
            fq2,
            star_align_file)

        cmd_bamindex_fltr = "{0} index {1}fltr_Aligned.sortedByCoord.out.bam".format(
            cmds["samtools"],
            star_align_file)

        cmd_requantify_fltr = "easy-fuse requantify " \
                              "-i {0}fltr_Aligned.sortedByCoord.out.bam " \
                              "-o {1}_fltr.tdt " \
                              "-d 10".format(
            star_align_file,
            classification_file)

        cmd_staralign_org = "{0} --genomeDir {1} --readFilesCommand zcat --readFilesIn {2} {3} " \
                            "--outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax -1 " \
                            "--outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 " \
                            "--runThreadN waiting_for_cpu_number --outFileNamePrefix {4}org_ " \
                            "--limitBAMsortRAM 48000000000".format(
            cmds["star"],
            star_genome_path,
            fq1,
            fq2,
            star_align_file)

        cmd_bamindex_org = "{0} index {1}org_Aligned.sortedByCoord.out.bam".format(
            cmds["samtools"],
            star_align_file)

        cmd_requantify_org = "easy-fuse requantify " \
                             "-i {0}org_Aligned.sortedByCoord.out.bam " \
                             "-o {1}_org.tdt -d 10".format(
            star_align_file,
            classification_file)

        # for testing, based on debug. should be removed if merged to original
        cmd_read_filter2 = "easy-fuse read-selection " \
                           "--input {0}_Aligned.out.bam " \
                           "--input2 {1}.debug " \
                           "--output {0}_Aligned.out.filtered2.bam".format(
            os.path.join(filtered_reads_path, sample_id),
            context_seq_file)

        # re-define fastq's if filtering is on (default)
        fq0 = ""
        if "readfilter" in tools:
            fq0 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(".fastq.gz", "_filtered2_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered2.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered2.fastq.gz"))
        cmd_bam_to_fastq_fd = "{0} fastq -0 {1} -1 {2} -2 {3} --threads waiting_for_cpu_number {4}_Aligned.out.filtered2.bam".format(cmds["samtools"], fq0, fq1, fq2, os.path.join(filtered_reads_path, sample_id))

        # allow soft-clipping? Specificity? --alignEndsType EndToEnd
        cmd_staralign_best = "{0} " \
                             "--genomeDir {1} " \
                             "--readFilesCommand zcat " \
                             "--readFilesIn {2} {3} " \
                             "--outSAMtype BAM SortedByCoordinate " \
                             "--outFilterMultimapNmax -1 " \
                             "--outSAMattributes Standard " \
                             "--outSAMunmapped None " \
                             "--outFilterMismatchNoverLmax 0.02 " \
                             "--runThreadN waiting_for_cpu_number " \
                             "--outFileNamePrefix {4}best_ " \
                             "--limitBAMsortRAM 48000000000".format(
            cmds["star"],
            star_genome_path,
            fq1,
            fq2,
            star_align_file)

        cmd_bamindex_best = "{0} index " \
                            "{1}best_Aligned.sortedByCoord.out.bam".format(
            cmds["samtools"],
            star_align_file)

        cmd_requantify_best = "easy-fuse requantify " \
                              "-i {0}best_Aligned.sortedByCoord.out.bam " \
                              "-o {1}_best.tdt -d 10".format(
            star_align_file,
            classification_file)

        # set final lists of executable tools and path
        exe_tools = [
            "qc", #0
            READ_FILTER_STEP, #1
            "star",
            "mapsplice",
            "fusioncatcher",
            "starfusion",
            "infusion",
            "soapfuse",
            "readcounts",
            "fusiongrep",
            "contextseq",
            "starindex",
            "staralignFltr",
            "bamindexFltr",
            "requantifyFltr",
            "staralignOrg",
            "bamindexOrg",
            "requantifyOrg",
            "readFilter2",
            "readFilter2b",
            "staralignBest",
            "bamindexBest",
            "requantifyBest"
            ]
        exe_cmds = [
            " && ".join([cmd_fastqc, cmd_qc_parser, cmd_skewer]),
            " && ".join([cmd_star_filter, cmd_read_filter, cmd_bam_to_fastq]),
            cmd_star,
            " && ".join([cmd_extr_fastq1, cmd_extr_fastq2, cmd_mapsplice]),
            cmd_fusioncatcher,
            cmd_starfusion,
            cmd_infusion,
            cmd_soapfuse,
            cmd_readcounts,
            cmd_fusiondata,
            cmd_contextseq,
            cmd_starindex,
            cmd_staralign_fltr,
            cmd_bamindex_fltr,
            cmd_requantify_fltr,
            cmd_staralign_org,
            cmd_bamindex_org,
            cmd_requantify_org,
            cmd_read_filter2,
            cmd_bam_to_fastq_fd,
            cmd_staralign_best,
            cmd_bamindex_best,
            cmd_requantify_best
            ]
        exe_path = [
            qc_path,
            filtered_reads_path,
            star_path,
            mapsplice_path,
            fusioncatcher_path,
            starfusion_path,
            infusion_path,
            soapfuse_path,
            classification_path,
            detected_fusions_path,
            context_seq_path,
            star_genome_path,
            star_align_path,
            star_align_path,
            classification_path,
            star_align_path,
            star_align_path,
            classification_path,
            filtered_reads_path,
            filtered_reads_path,
            star_align_path,
            star_align_path,
            classification_path
            ]

        # create and submit slurm job if the tool is requested and hasn't been run before
        for i, tool in enumerate(exe_tools, 0):
            if tool in tools:
                # check dependencies of the pipeline.
                # Besides tool dependencies (Pizzly -> Kallisto, Starfusion/Starchip -> Star), read filtering is mandatory
                if tool == READ_FILTER_STEP and READ_FILTER_STEP not in tools:
                    logger.error(
                            """Error 99: Sample {} will be skipped due to missing read filtering.\n
                            Read filtering is currently a mandatory step for the processing.\n
                            Because you haven't run it before for this sample, you have to include \"Readfilter\" in the tool selection in your config.\n
                            """.format(sample_id))
                    print("Error 99: Sample {} will be skipped due to missing read filtering.".format(sample_id))
                    return 0
                elif tool == "pizzly" and "kallisto" not in tools:
                    logger.error(
                            """Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.\n
                            Pizzly builds on Kallisto and it is therefore mandatory to run this first.\n
                            Because you haven't run it before for this sample, you have to include \"Kallisto\" in the tool selection in your config.\n
                            """.format(sample_id))
                    print("Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.".format(tool, sample_id))
                    continue
                elif (tool == "starfusion" or tool == "starchip") and "star" not in tools:
                    logger.error(
                            """Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.\n
                            {0} builds on Star and it is therefore mandatory to run this first.\n
                            Because you haven't run it before for this sample, you have to include \"Star\" in the tool selection in your config.\n
                            """.format(sample_id))
                    print("Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.".format(tool, sample_id))
                    continue

                # prepare slurm jobs: get ressources, create uid, set output path and check dependencies
                logger.debug("Submitting {} run".format(tool))
                cpumem = self.cfg["resources"].get(tool.lower(), DEFAULT_CPU_MEM).split(",")
                cpu = cpumem[0]
                mem = cpumem[1]
                uid = "-".join([self.cfg["general"]["pipeline_name"], tool, sample_id])
                if tool == "star":
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string", "{}_".format(os.path.join(exe_path[i], sample_id))).replace("waiting_for_cpu_number", str(cpu))
                else:
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string", exe_path[i]).replace(
                        "waiting_for_cpu_number", str(cpu))
                cmd = exe_cmds[i]

                # Build slurm dependencies
                dependencies = self.build_dependencies(tool, sample_id)

                logger.debug("Submitting job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd, exe_path[i], dependencies))
                self.submit_job(uid, cmd, cpu, mem, exe_path[i], dependencies, "")
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



    def build_dependencies(self, tool, sample_id):
        slurm_dep = {
            "qc": [],
            "readfilter": ["qc"],
            "star": ["readfilter"],
            "readcounts": ["star"],
            "starfusion": ["star"],
            "mapsplice": ["readfilter"],
            "fusioncatcher": ["readfilter"],
            "infusion": ["readfilter"],
            "soapfuse": ["readfilter"],
            "fusiongrep": ["starfusion", "mapsplice", "fusioncatcher", "infusion", "soapfuse"],
            "contextseq": ["fusiongrep"],
            "starindex": ["readcounts", "contextseq"],
            "readFilter2": ["contextseq"],
            "readFilter2b": ["readFilter2"],
            "staralignBest": ["starindex", "readFilter2b"],
            "bamindexBest": ["staralignBest"],
            "requantifyBest": ["bamindexBest"]
        }
        dependencies = []
        que_sys = self.cfg["general"]["queueing_system"]

        for pre_tool in slurm_dep[tool]:
            dependencies.extend(queueing.get_jobs_by_name("{}-{}".format(pre_tool, sample_id), que_sys))

        return dependencies


def add_pipeline_parser_args(pipeline_parser):
    pipeline_parser.add_argument(
        "-i",
        "--input",
        dest="input_paths",
        nargs="+",
        help="Specify full path of the fastq folder to process.",
        required=True,
    )
    pipeline_parser.add_argument(
        "-o",
        "--output-folder",
        dest="output_folder",
        help="Specify full path of the folder to save the results into.",
        required=True,
    )
    # optional arguments
    pipeline_parser.add_argument(
        "-c",
        "--config-file",
        dest="config_file",
        required=True,
        help="Specify alternative config file to use for your analysis",
    )
    pipeline_parser.add_argument(
        "-p",
        "--jobname_suffix",
        dest="jobname_suffix",
        help="Specify a jobname suffix for running jobs on a queueing system",
        default="",
    )
    pipeline_parser.add_argument(
        "-V",
        "--version",
        dest="version",
        action="version",
        version="%(prog)s {version}".format(version=easy_fuse.__version__),
        help="Get version info",
    )
    pipeline_parser.set_defaults(func=pipeline_command)


def pipeline_command(args):
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
        os.path.abspath(args.output_folder),
    )
    with open(os.path.join(args.output_folder, "process.sh"), "w") as outf:
        outf.write("#!/bin/sh\n\n{}".format(script_call))

    Processing(script_call, args.input_paths, args.output_folder, config, args.jobname_suffix).run()
