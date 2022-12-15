#!/usr/bin/env python3

"""
Read fastq files
create per-sample (sub)folder tree
run gene fusion predictions
start data collection
during all steps, perform slurm scheduling and process monitoring

@author: Tron (PASO), BNT (URLA)
"""


from argparse import ArgumentParser
from configparser import ConfigParser
import json
import os
import os.path
from shutil import copy
import sys
import time


from _version import __version__
import misc.queueing as Queueing
from misc.samples import SamplesDB
from misc.logger import Logger
import misc.io_methods as IOMethods


class Processing(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""
    def __init__(self, cmd, input_paths, working_dir, cfg_file, jobname_suffix, que_sys, partition):
        """Parameter initiation and work folder creation. Start of progress logging."""
        self.working_dir = os.path.abspath(working_dir)
        self.logger = Logger(os.path.join(self.working_dir, "easyfuse_processing.log"))
        IOMethods.create_folder(self.working_dir, self.logger)

        self.logger.info("Starting easyfuse: CMD - {}".format(cmd))
        self.input_paths = [os.path.abspath(file) for file in input_paths]
        self.samples = SamplesDB(os.path.join(self.working_dir, "samples.db"))

        self.cfg = None

        default_cfg_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini")

        if not cfg_file and os.path.exists(default_cfg_path):
            cfg_file = default_cfg_path
        elif not cfg_file and not os.path.exists(default_cfg_path):
            self.logger.error("Could not find default config file (path={}). Exiting.".format(default_cfg_path))
            sys.exit(1)

        self.cfg_file = cfg_file

        if cfg_file.endswith("ini"):
            self.cfg = ConfigParser()
            self.cfg.read(cfg_file)
        elif cfg_file.endswith("json"):
            with open(cfg_file) as config_file:
                self.cfg = json.load(config_file)

        copy(cfg_file, working_dir)


        self.jobname_suffix = None
        if jobname_suffix:
            self.jobname_suffix = jobname_suffix

        self.que_sys = que_sys
        self.partition = partition

    # The run method simply greps and organises fastq input files.
    # Fastq pairs (single end input is currently not supported) are then send to "execute_pipeline"
    def run(self):
        """General parameter setting, identification of fastq files and initiation of processing"""
        self.logger.info("Pipeline Version: {}".format(__version__))

        ref_genome = self.cfg["general"]["ref_genome_build"]
        ref_trans = self.cfg["general"]["ref_trans_version"]

        self.logger.info("Reference Genome: {0}, Reference Transcriptome: {1}".format(ref_genome, ref_trans))


        # get fastq files
        fastqs = IOMethods.get_fastq_files(self.input_paths, self.logger)
        sample_list = IOMethods.pair_fastq_files(fastqs, self.logger)
        for sample in sample_list:
            if sample[1] and sample[2]:
                self.logger.info("Processing Sample ID: {} (paired end)".format(sample[0]))
                self.logger.info("Sample 1: {}".format(sample[1]))
                self.logger.info("Sample 2: {}".format(sample[2]))
                self.execute_pipeline(sample[1], sample[2], sample[0], ref_genome, ref_trans)

        
        # summarize all data if selected
        if "summary" in self.cfg["general"]["tools"].split(","):
            dependency = []
            for sample in sample_list:
                dependency.extend(Queueing.get_jobs_by_name("requantifyBest-{}".format(sample[0]), self.que_sys))
            modelling_string = ""
            if self.cfg["other_files"]["easyfuse_model"]:
                modelling_string = " --model_predictions"
            cmd_summarize = "{0} {1} --input {2}{3} -c {4}".format(self.cfg["commands"]["python3"], os.path.join(self.cfg["general"]["module_dir"], "summarize_data.py"), self.working_dir, modelling_string, self.cfg_file)
            self.logger.debug("Submitting job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd_summarize, self.working_dir, dependency))
            resources_summary = self.cfg["resources"]["summary"].split(",")
            cpu = resources_summary[0]
            mem = resources_summary[1]
            self.submit_job("-".join([self.cfg["general"]["pipeline_name"], "Summary", str(int(round(time.time())))]), cmd_summarize, cpu, mem, self.working_dir, dependency, self.cfg["general"]["receiver"])

    # Per sample, define input parameters and execution commands, create a folder tree and submit runs to slurm
    def execute_pipeline(self, fq1, fq2, sample_id, ref_genome, ref_trans):
        """Create sample specific subfolder structure and run tools on fastq files"""
        self.samples.add_sample(sample_id, "NA", fq1, fq2)
                
        refs = self.cfg["references"]

        # Genome/Gene references to use
        genome_sizes_path = refs["genome_sizes"]
        genome_chrs_path = refs["genome_fastadir"]
        genome_fasta_path = refs["genome_fasta"]
        
        genes_fasta_path = refs["genes_fasta"]
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
        fastqc_1 = os.path.join(qc_path, os.path.basename(fq1).rstrip(".fastq.gz") + "_fastqc", "fastqc_data.txt")
        fastqc_2 = os.path.join(qc_path, os.path.basename(fq2).rstrip(".fastq.gz") + "_fastqc", "fastqc_data.txt")
        skewer_path = os.path.join(qc_path, "skewer")
        qc_table_path = os.path.join(qc_path, "qc_table.txt")
        overrepresented_path = os.path.join(qc_path, "overrepresented.fa")
        filtered_reads_path = os.path.join(output_results_path, "filtered_reads")
        expression_path = os.path.join(output_results_path, "expression")
        star_path = os.path.join(expression_path, "star")
        fusion_path = os.path.join(output_results_path, "fusion")
        mapsplice_path = os.path.join(fusion_path, "mapsplice")
        pizzly_path = os.path.join(fusion_path, "pizzly")
        fusioncatcher_path = os.path.join(fusion_path, "fusioncatcher")
        starfusion_path = os.path.join(fusion_path, "starfusion")
        starchip_path = os.path.join(fusion_path, "starchip")
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
            IOMethods.create_folder(folder, self.logger)

        # get a list of tools from the samples.db file that have been run previously on this sample
        state_tools = self.samples.get_tool_list_from_state(sample_id)
        #print(sample_id, state_tools)
        # get a list of tools from the config file which shall be run on this sample
        tools = self.cfg["general"]["tools"].split(",")
        cmds = self.cfg["commands"]
        module_dir = self.cfg["general"]["module_dir"]
        # Define cmd strings for each program
        # (0) QC
        cmd_fastqc = "{0} --nogroup --extract -t 6 -o {1} {2} {3}".format(cmds["fastqc"], qc_path, fq1, fq2)
        cmd_qc_parser = "{0} {1} -i {2} {3} -o {4}".format(cmds["python3"], os.path.join(module_dir, "misc", "qc_parser.py"), fastqc_1, fastqc_2, qc_table_path)
        cmd_skewer = "{0} {1} -q {2} -i {3} {4} -o {5} -b {6} -m {7}".format(cmds["python3"], os.path.join(module_dir, "tool_wrapper", "skewer_wrapper.py"), qc_table_path, fq1, fq2, skewer_path, self.cfg["commands"]["skewer"], self.cfg["general"]["min_read_len_perc"])

        fq0 = ""
        if "qc" in tools:
            fq0 = os.path.join(skewer_path, "out_file-trimmed.fastq.gz")
            fq1 = os.path.join(skewer_path, "out_file-trimmed-pair1.fastq.gz")
            fq2 = os.path.join(skewer_path, "out_file-trimmed-pair2.fastq.gz")
        else:
            qc_table_path = "None"

        # (1) Readfilter
        cmd_star_filter = "{0} --genomeDir {1} --outFileNamePrefix {2}_ --readFilesCommand zcat --readFilesIn {3} {4} --outFilterMultimapNmax 100 --outSAMmultNmax 1 --chimSegmentMin 10 --chimJunctionOverhangMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax {5} --alignIntronMax {5} --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --seedSearchStartLmax 20 --winAnchorMultimapNmax 50 --outSAMtype BAM Unsorted --chimOutType Junctions WithinBAM --outSAMunmapped Within KeepPairs --runThreadN waiting_for_cpu_number".format(cmds["star"], star_index_path, os.path.join(filtered_reads_path, sample_id), fq1, fq2, self.cfg["general"]["max_dist_proper_pair"])
        cmd_read_filter = "{0} {1} --input {2}_Aligned.out.bam --output {2}_Aligned.out.filtered.bam".format(cmds["python3"], os.path.join(module_dir, "fusionreadfilter.py"), os.path.join(filtered_reads_path, sample_id))
        # re-define fastq's if filtering is on (default)
        fq0 = ""
        if "readfilter" in tools:
            fq0 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(".fastq.gz", "_filtered_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered.fastq.gz"))

        cmd_bam_to_fastq = "{0} fastq -0 {1} -1 {2} -2 {3} --threads waiting_for_cpu_number {4}_Aligned.out.filtered.bam".format(cmds["samtools"], fq0, fq1, fq2, os.path.join(filtered_reads_path, sample_id))
        # (2) Star expression quantification (required for starfusion)
        cmd_star = "{0} --genomeDir {1} --outFileNamePrefix waiting_for_output_string --runThreadN waiting_for_cpu_number --runMode alignReads --readFilesIn {2} {3} --readFilesCommand zcat --chimSegmentMin 10 --chimJunctionOverhangMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax {4} --alignIntronMax {4} --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --seedSearchStartLmax 20 --winAnchorMultimapNmax 50 --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold --chimOutJunctionFormat 1".format(cmds["star"], star_index_path, fq1, fq2, self.cfg["general"]["max_dist_proper_pair"])
        # (3) Mapslice
        # Using "gunzip -c" instead of "gunzip --keep" to keep compatibility with gunzip versions
        cmd_extr_fastq1 = "gunzip -c -f {0} > {1}".format(fq1, fq1[:-3])
        cmd_extr_fastq2 = "gunzip -c -f {0} > {1}".format(fq2, fq2[:-3])
        # Added python interpreter to circumvent external hardcoded shell script
        cmd_mapsplice = "{0} --chromosome-dir {1} -x {2} -1 {3} -2 {4} --threads waiting_for_cpu_number --output {5} --qual-scale phred33 --bam --seglen 20 --min-map-len 40 --gene-gtf {6} --fusion".format(cmds["mapsplice"], genome_chrs_path, bowtie_index_path, fq1[:-3], fq2[:-3], mapsplice_path, genes_gtf_path)
        # (4) Fusioncatcher
        cmd_fusioncatcher = "{0} --input {1} --data {2} --output {3} -p waiting_for_cpu_number".format(cmds["fusioncatcher"], ",".join([fq1, fq2]), fusioncatcher_index_path, fusioncatcher_path)
        # (5) Starfusion
        cmd_starfusion = "{0} --chimeric_junction {1} --genome_lib_dir {2} --CPU waiting_for_cpu_number --output_dir {3}".format(cmds["starfusion"], os.path.join(star_path, "{}_Chimeric.out.junction".format(sample_id)), starfusion_index_path, starfusion_path)
        # (6) Infusion
        # TODO: Use these command line args for more sensitive results:
        # --min-split-reads 0 --min-fragments 0 --min-not-rescued 2
        cmd_infusion = "{0} {1} -1 {2} -2 {3} --skip-finished --min-unique-alignment-rate 0 --min-unique-split-reads 0 --allow-non-coding --out-dir {4} {5}".format(cmds["python2"], cmds["infusion"], fq1, fq2, infusion_path, infusion_cfg_path)
        # (7) Soapfuse
        cmd_soapfuse = "{0} {1} -q {2} -i {3} -o {4} -b {5} -c {6}".format(cmds["python3"], os.path.join(module_dir, "tool_wrapper", "soapfuse_wrapper.py"), qc_table_path, " ".join([fq1, fq2]), soapfuse_path, self.cfg["commands"]["soapfuse"], self.cfg["other_files"]["soapfuse_cfg"])
        # (8) Data collection
        cmd_readcounts = "{0} {1} -i {2} -f {3} -o {4}".format(cmds["python3"], os.path.join(module_dir, "misc", "count_input_reads.py"), os.path.join(filtered_reads_path, "{}_Log.final.out".format(sample_id)), "star", os.path.join(classification_path, "Star_org_input_reads.txt"))
        cmd_fusiondata = "{0} {1} -i {2} -o {3} -s {4} -f {5} -l {6}".format(cmds["python3"], os.path.join(module_dir, "fusiontoolparser.py"), output_results_path, detected_fusions_path, sample_id, self.cfg["general"]["fusiontools"], self.logger.get_path())
        cmd_contextseq = "{0} {1} --detected_fusions {2} --annotation_db {3} --out_csv {4} --genome_fasta {5} --tsl_info {6} --cis_near_dist {7} --context_seq_len {8} --tsl_filter_level {9}".format(cmds["python3"], os.path.join(module_dir, "fusionannotation.py"), detected_fusions_file, genes_adb_path, context_seq_file, genome_fasta_path, genes_tsl_path, self.cfg["general"]["cis_near_distance"], self.cfg["general"]["context_seq_len"], self.cfg["general"]["tsl_filter"])
        cmd_starindex = "{0} {1} -i {2} -o {3} -t waiting_for_cpu_number -b {4}".format(cmds["python3"], os.path.join(module_dir, "tool_wrapper", "star_custom_index.py"), context_seq_fasta, star_genome_path, cmds["star"])
        cmd_staralign_fltr = "{0} --genomeDir {1} --readFilesCommand zcat --readFilesIn {2} {3} --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax -1 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 --runThreadN waiting_for_cpu_number --outFileNamePrefix {4}fltr_ --limitBAMsortRAM 48000000000".format(cmds["star"], star_genome_path, fq1, fq2, star_align_file)
        cmd_bamindex_fltr = "{0} index {1}fltr_Aligned.sortedByCoord.out.bam".format(cmds["samtools"], star_align_file)
        cmd_requantify_fltr = "{0} -i {1}fltr_Aligned.sortedByCoord.out.bam -o {2}_fltr.tdt -d 10".format(os.path.join(module_dir, "requantify.py"), star_align_file, classification_file)
        (fq1, fq2) = self.samples.get_fastq_files(sample_id)
        cmd_staralign_org = "{0} --genomeDir {1} --readFilesCommand zcat --readFilesIn {2} {3} --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax -1 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 --runThreadN waiting_for_cpu_number --outFileNamePrefix {4}org_ --limitBAMsortRAM 48000000000".format(cmds["star"], star_genome_path, fq1, fq2, star_align_file)
        cmd_bamindex_org = "{0} index {1}org_Aligned.sortedByCoord.out.bam".format(cmds["samtools"], star_align_file)
        cmd_requantify_org = "{0} {1} -i {2}org_Aligned.sortedByCoord.out.bam -o {3}_org.tdt -d 10".format(cmds["python3"], os.path.join(module_dir, "requantify.py"), star_align_file, classification_file)
        # for testing, based on debug. should be removed if merged to original
        cmd_read_filter2 = "{0} {1} --input {2}_Aligned.out.bam --input2 {3}.debug --output {2}_Aligned.out.filtered2.bam".format(cmds["python3"], os.path.join(module_dir, "getRequantReads.py"), os.path.join(filtered_reads_path, sample_id), context_seq_file)

        # re-define fastq's if filtering is on (default)
        fq0 = ""
        if "readfilter" in tools:
            fq0 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(".fastq.gz", "_filtered2_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered2.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered2.fastq.gz"))
        cmd_bam_to_fastq_fd = "{0} fastq -0 {1} -1 {2} -2 {3} --threads waiting_for_cpu_number {4}_Aligned.out.filtered2.bam".format(cmds["samtools"], fq0, fq1, fq2, os.path.join(filtered_reads_path, sample_id))
        # allow soft-clipping? Specificity? --alignEndsType EndToEnd
        cmd_staralign_best = "{0} --genomeDir {1} --readFilesCommand zcat --readFilesIn {2} {3} --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax -1 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 --runThreadN waiting_for_cpu_number --outFileNamePrefix {4}best_ --limitBAMsortRAM 48000000000".format(cmds["star"], star_genome_path, fq1, fq2, star_align_file)
        cmd_bamindex_best = "{0} index {1}best_Aligned.sortedByCoord.out.bam".format(cmds["samtools"], star_align_file)
        cmd_requantify_best = "{0} {1} -i {2}best_Aligned.sortedByCoord.out.bam -o {3}_best.tdt -d 10".format(cmds["python3"], os.path.join(module_dir, "requantify.py"), star_align_file, classification_file)
        
        # (X) Sample monitoring
        cmd_samples = "{0} {1} --db_path={2} --sample_id={3} --action=append_state --tool=".format(cmds["python3"], os.path.join(module_dir, "misc", "samples.py"), self.samples.db_path, sample_id)

        # set final lists of executable tools and path
        exe_tools = [
            "qc",
            "readfilter",
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
                dependency = []
                # check dependencies of the pipeline.
                # Besides tool dependencies (Pizzly -> Kallisto, Starfusion/Starchip -> Star), read filtering is mandatory
                # Processing will be skipped if a certain dependency was not found (either pre-processed data of the configs tool string are checked)
                if tool in state_tools:
                    self.logger.info("Skipping {0} as it looks like a previous run finished successfully. Results should be in {1}".format(tool, exe_path[i]))
                    continue
                else:
                    if tool == "readfilter" and "readfilter" not in tools:
                        self.logger.error(
                                """Error 99: Sample {} will be skipped due to missing read filtering.\n
                                Read filtering is currently a mandatory step for the processing.\n
                                Because you haven't run it before for this sample, you have to include \"Readfilter\" in the tool selection in your config.\n
                                """.format(sample_id))
                        print("Error 99: Sample {} will be skipped due to missing read filtering.".format(sample_id))
                        return 0
                    elif tool == "pizzly" and "kallisto" not in tools:
                        self.logger.error(
                                """Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.\n
                                Pizzly builds on Kallisto and it is therefore mandatory to run this first.\n
                                Because you haven't run it before for this sample, you have to include \"Kallisto\" in the tool selection in your config.\n
                                """.format(sample_id))
                        print("Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.".format(tool, sample_id))
                        continue
                    elif (tool == "starfusion" or tool == "starchip") and "star" not in tools:
                        self.logger.error(
                                """Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.\n
                                {0} builds on Star and it is therefore mandatory to run this first.\n
                                Because you haven't run it before for this sample, you have to include \"Star\" in the tool selection in your config.\n
                                """.format(sample_id))
                        print("Error 99: Running {0} for Sample {1} will be skipped due to a missing dependency.".format(tool, sample_id))
                        continue

                # prepare slurm jobs: get ressources, create uid, set output path and check dependencies
                self.logger.debug("Submitting {} run".format(tool))
                cpumem = self.cfg["resources"][tool.lower()].split(",")
                cpu = cpumem[0]
                mem = cpumem[1]
                uid = "-".join([self.cfg["general"]["pipeline_name"], tool, sample_id])
                if tool == "star":
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string", "{}_".format(os.path.join(exe_path[i], sample_id))).replace("waiting_for_cpu_number", str(cpu))
                else:
                    exe_cmds[i] = exe_cmds[i].replace("waiting_for_output_string", exe_path[i]).replace("waiting_for_cpu_number", str(cpu))
                cmd = " && ".join([exe_cmds[i], cmd_samples + tool])
                
                # Build slurm dependencies
                dependencies = self.build_dependencies(tool, sample_id)

                self.logger.debug("Submitting job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd, exe_path[i], dependencies))
                self.submit_job(uid, cmd, cpu, mem, exe_path[i], dependencies, "")
            else:
                self.logger.info("Skipping {0} as it is not selected for execution (Selected are: {1})".format(tool, tools))


    def submit_job(self, uid, cmd, cores, mem_usage, output_results_folder, dependencies, mail):
        """Submit job to for process generation"""
        already_running = Queueing.get_jobs_by_name(uid, self.que_sys)
        if not already_running:
            module_file = os.path.join(self.cfg["general"]["module_dir"], "build_env.sh")
            # for compatibility reasons (and to be independent of shell commands), 
            # concatenated commands are splitted again,
            # dependencies within the splitted groups updated 
            # and everything submitted sequentially to the queueing system
            for i, cmd_split in enumerate(cmd.split(" && ")):
                if not self.que_sys in ["slurm", "pbs"]:
                    cmd_split = cmd_split.split(" ")
                dependencies.extend(Queueing.get_jobs_by_name("{0}_CMD{1}".format(uid, i - 1), self.que_sys))
                if self.jobname_suffix:
                    Queueing.submit("{0}_CMD{1}-{2}".format(uid, i, self.jobname_suffix), cmd_split, cores, mem_usage, output_results_folder, dependencies, self.partition, os.getlogin(), self.cfg["general"]["time_limit"], mail, module_file, self.que_sys)
                else:
                    Queueing.submit("{0}_CMD{1}".format(uid, i), cmd_split, cores, mem_usage, output_results_folder, dependencies, self.partition, os.getlogin(), self.cfg["general"]["time_limit"], mail, module_file, self.que_sys)
                time.sleep(0.5)
        else:
            self.logger.error("A job with this application/sample combination is currently running. Skipping {} in order to avoid unintended data loss.".format(uid))


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

        for pre_tool in slurm_dep[tool]:
            dependencies.extend(Queueing.get_jobs_by_name("{}-{}".format(pre_tool, sample_id), self.que_sys))

        return dependencies


def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(prog='EasyFuse', description='Run EasyFuse pipeline')
    # required arguments
    parser.add_argument(
        '-i',
        '--input',
        dest='input_paths',
        nargs='+',
        help='Specify full path of the fastq folder to process.',
        required=True
    )
    parser.add_argument(
        '-o',
        '--output_folder',
        dest='output_folder',
        help='Specify full path of the folder to save the results into.',
        required=True
    )
    # optional arguments
    parser.add_argument(
        '-c',
        '--config_file',
        dest='config_file',
        help='Specify alternative config file to use for your analysis'
    )
    parser.add_argument(
        '-p',
        '--partition',
        dest='partition',
        help='Specify partition for job submission in slurm, where partitions can be viewed using "sinfo"'
    )
    parser.add_argument(
        '-q',
        '--queueing_system',
        dest='queueing_system',
        choices=['none', 'slurm'],
        help='Specify queueing system to be used',
        default='none'
    )
    parser.add_argument(
        '-s',
        '--jobname_suffix',
        dest='jobname_suffix',
        help='Specify a jobname suffix for running jobs on a queueing system',
        default=''
    )
    parser.add_argument(
        '-V',
        '--version',
        dest='version',
        action='version',
        version="%(prog)s {version}".format(version=__version__),
        help='Get version info'
    )
    args = parser.parse_args()

    jobname_suffix = ""
    if args.jobname_suffix:
        jobname_suffix = "-p {}".format(args.jobname_suffix)

    config = ""
    cfg_file = None
    if args.config_file:
        cfg_file = os.path.abspath(args.config_file)
        config = "-c {}".format(cfg_file)

    script_call = "python {} -i {} {} {} -o {}".format(
        os.path.realpath(__file__),
        " ".join([os.path.abspath(x) for x in args.input_paths]),
        config,
        jobname_suffix,
        os.path.abspath(args.output_folder)
    )

    proc = Processing(script_call, args.input_paths, args.output_folder, cfg_file, args.jobname_suffix, args.queueing_system, args.partition)
    proc.run()

    with open(os.path.join(args.output_folder, "process.sh"), "w") as outf:
        outf.write("#!/bin/sh\n\n{}".format(script_call))

if __name__ == '__main__':
    main()
