"""
Read fastq files
create per-sample (sub)folder tree
run gene fusion predictions
start data collection
during all steps, perform slurm scheduling and process monitoring

@author: Tron (PASO), BNT (URLA)
@version: 20190118
"""
from __future__ import print_function
import os
import os.path
import sys
import time
import argparse
from shutil import copy2
import misc.queue as Queueing
from misc.config import Config
from misc.samples import Samples
from misc.logger import Logger
import misc.io_methods as IOMethods
from misc.versioncontrol import VersCont
import misc.version as eav

# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class Processing(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""
    def __init__(self, config, input_path, working_dir, partitions, userid, overwrite):
        """Parameter initiation and work folder creation. Start of progress logging."""
        self.working_dir = working_dir.rstrip("/")
        self.logger = Logger(os.path.join(self.working_dir, "easyfuse_processing_{}.log".format(str(int(round(time.time()))))))
        IOMethods.create_folder(self.working_dir, self.logger)
        copy2(config, working_dir)
        self.cfg = Config(os.path.join(self.working_dir, IOMethods.path_leaf(config)))
        self.logger.info("Starting easyfuse")
        self.input_path = input_path
        self.samples = Samples(os.path.join(self.working_dir, "samples.csv"))
        self.partitions = partitions
        self.userid = userid
        self.overwrite = overwrite

    # The run method simply greps and organises fastq input files.
    # Fastq pairs (single end input is currently not supported) are then send to "execute_pipeline"
    def run(self, tool_num_cutoff, run_qc, filter_reads, icam_run, icam_tree):
        """General parameter setting, identification of fastq files and initiation of processing"""
        # Checking dependencies
        VersCont(self.cfg.get('easyfuse_helper', 'dependencies')).get_and_print_tool_versions()
        self.cfg.run_self_test()
        is_stranded = False # urla: I leave it in although it is currently not required - however, future predictions may improve from stranded libs
        is_rna_seq = True # urla: can this ever become false?!
        # urla: organism is currently not used for anything, however, this might change; is mouse processing relevant at some point?
        ref_genome = self.cfg.get('general', 'ref_genome_build')
        ref_trans = self.cfg.get('general', 'ref_trans_version')
#        if ref_genome in ['hg19', 'hg38']:
#            organism = "human"
#        else:
#            organism = "mouse"
        self.logger.info("Stranded: {0}, RNA-Seq: {1}, Reference Genome: {2}, Reference Transcriptome: {3}".format(is_stranded, is_rna_seq, ref_genome, ref_trans))
        if self.overwrite:
            self.logger.info("#############################################################################")
            self.logger.info("")
            self.logger.info("Overwrite flag is set => all previously existing results may be overwritten!")
            self.logger.info("")
            self.logger.info("#############################################################################")

        # get fastq files
        # urla: the icam_run flag starts a file walker which searches fastq files in a dir-tree containing icam-flowcell folders
        #       An example of an input root folder would be "/mnt/bfx/IVAC2_0/BNT/iCaM2Scheduler/data/demux_out"
        #       Only folders containing a valid (see iomethods for details) "demux_stats.csv" file can be processed atm
        # urla: the icam_run flag is probably not necessary for production usage but very convenient if novel tools/methods should be run on a bunch of existing datasets (future development / evaluations)
        sample_list = []
        if icam_run and icam_tree:
            self.logger.info("Started processing of icam data; searching fastq files in \"{}\"".format(self.input_path))
            # print(len(self.samples.sample_map)) <- 0 for empty list
            _, fastq_dict, counter_list = IOMethods.get_icam_reads_data(self.input_path, list(self.samples.sample_map), self.logger)
            count_all = sum(counter_list)
            self.logger.debug("Found a total of {0} datasets in the specified folder, {1} have been processed previously, {2} are new.".format(count_all, counter_list[0], counter_list[1]))
            self.logger.debug("{0} datasets could not be processed due to an unexpected format of the demux_stats.csv file.".format(counter_list[2]))
            self.logger.debug("{0} datasets were not be processed because their read count was below 75M.".format(counter_list[3]))
            self.logger.debug("{0} datasets could not be processed due problems in the fastq file localization.".format(counter_list[4]))
            for i, sample_id_key in enumerate(sorted(fastq_dict)):
                sample_id = "_".join(sample_id_key.split("_")[:2])
                sample_list.append(sample_id)
                self.logger.info("Processing Sample ID: {} (paired end)".format(sample_id))
                self.logger.info("Sample 1: {}".format(fastq_dict[sample_id_key][0]))
                self.logger.info("Sample 2: {}".format(fastq_dict[sample_id_key][1]))
                self.execute_pipeline(fastq_dict[sample_id_key][0], fastq_dict[sample_id_key][1], sample_id, ref_genome, ref_trans, tool_num_cutoff, run_qc, filter_reads, icam_run)
        else:
            left, right, sample_id = IOMethods.get_fastq_files(self.input_path, self.logger)
            sample_list = sample_id
            for i, _ in enumerate(left):
                if len(left) == len(right):
                    self.logger.info("Processing Sample ID: {} (paired end)".format(sample_id[i]))
                    self.logger.info("Sample 1: {}".format(left[i]))
                    self.logger.info("Sample 2: {}".format(right[i]))
                    self.execute_pipeline(left[i], right[i], sample_id[i], ref_genome, ref_trans, tool_num_cutoff, run_qc, filter_reads, icam_run)
        
        # summarize all data if selected
        if "Summary" in self.cfg.get('general', 'tools'):
            #dependency = [Queueing.get_jobs_by_name("Fetchdata-{}".format(sample)) for sample in sample_list]
            # urla - note: would be happy to get the dependencies with a stacked LC, but is atm to complicated for me ^^
            dependency = []
            for sample in sample_list:
                dependency.extend(Queueing.get_jobs_by_name("Fetchdata-{}".format(sample)))
            icam_run_string = ""
            if icam_run:
                icam_run_string = " --icam_run"
            modelling_string = ""
            if self.cfg.get("otherFiles", "easyfuse_model"):
                modelling_string = " --model_predictions"
            cmd_summarize = "{0} --input {1} --config {2}{3}{4}".format(self.cfg.get('commands', 'summary_cmd'), self.working_dir, self.cfg.get_path(), icam_run_string, modelling_string)
            self.logger.debug("Submitting slurm job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd_summarize, self.working_dir, dependency))
            cpu, mem = self.cfg.get("resources", "summary").split(",")
            self.submit_job("-".join(["Summary", str(int(round(time.time())))]), cmd_summarize, cpu, mem, self.working_dir, dependency)

    # Per sample, define input parameters and execution commands, create a folder tree and submit runs to slurm
    def execute_pipeline(self, fq1, fq2, sample_id, ref_genome, ref_trans, tool_num_cutoff, run_qc, filter_reads, icam_run):
        """Create sample specific subfolder structure and run tools on fastq files"""
        self.samples.add_sample(sample_id, "NA", "R1:{0};R2:{1}".format(fq1, fq2))

		# Genome/Gene references to use
        genome_sizes_path = self.cfg.get("references", "{0}_genome_sizes_{1}".format(ref_trans, ref_genome))
        genome_chrs_path = self.cfg.get("references", "{0}_genome_fastadir_{1}".format(ref_trans, ref_genome))
        genes_fasta_path = self.cfg.get("references", "{0}_genes_fasta_{1}".format(ref_trans, ref_genome))
        genes_gtf_path = self.cfg.get("references", "{0}_genes_gtf_{1}".format(ref_trans, ref_genome))

        # Path' to specific indices
        bowtie_index_path = self.cfg.get("indices", "{0}_bowtie_{1}".format(ref_trans, ref_genome))
        star_index_path = self.cfg.get("indices", "{0}_star_{1}_sjdb{2}".format(ref_trans, ref_genome, 49))
#        star_index_path = self.cfg.get("indices", "{0}_star_{1}_trans".format(ref_trans, ref_genome))
        kallisto_index_path = self.cfg.get("indices", "{0}_kallisto_{1}".format(ref_trans, ref_genome))
        pizzly_cache_path = "{}.pizzlyCache.txt".format(genes_gtf_path)
        starfusion_index_path = self.cfg.get("indices", "{0}_starfusion_{1}".format(ref_trans, ref_genome))
        infusion_cfg_path = self.cfg.get("otherFiles", "{0}_infusion_cfg_{1}".format(ref_trans, ref_genome))
        starchip_param_path = self.cfg.get("otherFiles", "{0}_starchip_param_{1}".format(ref_trans, ref_genome))

        # Output results folder creation - currently included:
        # 1) Gene/Isoform expression: kallisto, star
        # 2) Fusion prediction: mapsplice, pizzly, fusioncatcher, star-fusion, starchip, infusion
        output_results_path = os.path.join(self.working_dir, "Sample_{}".format(sample_id), "scratch")
        qc_path = os.path.join(output_results_path, "qc")
        skewer_path = os.path.join(qc_path, "skewer")
        qc_table_path = os.path.join(qc_path, "qc_table.txt")
        overrepresented_path = os.path.join(qc_path, "overrepresented.fa")
        filtered_reads_path = os.path.join(output_results_path, "filtered_reads")
        expression_path = os.path.join(output_results_path, "expression")
        kallisto_path = os.path.join(expression_path, "kallisto")
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

        for folder in [
                output_results_path, qc_path, skewer_path, filtered_reads_path,
                expression_path, kallisto_path, star_path,
                fusion_path, mapsplice_path, pizzly_path, fusioncatcher_path, starfusion_path, starchip_path, infusion_path, soapfuse_path,
                fetchdata_path
            ]:
            IOMethods.create_folder(folder, self.logger)

        # get a list of tools from the samples.csv file that have been run previously on this sample
        state_tools = self.samples.get_tool_list_from_state(sample_id)
        # get a list of tools from the config file which shall be run on this sample
        tools = []
        if "," in self.cfg.get('general', 'tools'):
            tools.extend(self.cfg.get('general', 'tools').split(","))
        else:
            tools.extend([self.cfg.get('general', 'tools')])

        # Define cmd strings for each program
        # urla: mapsplice requires gunzip'd read files and process substitutions don't seem to work in slurm scripts...
        #       process substitution do somehow not work from this script - c/p the command line to the terminal, however, works w/o issues?!

        cmd_fastqc = "{} --nogroup --extract -t 6 -o {} {} {}".format(self.cfg.get('commands', 'fastqc_cmd'), qc_path, fq1, fq2)
        cmd_qc_parser = "{} -i {}/*/fastqc_data.txt -q {} -o {}".format(self.cfg.get('commands', 'qc_parser_cmd'), qc_path, qc_table_path, overrepresented_path)
        cmd_skewer = "{} -q {} -m {} -i {} {} -o {} -b {}".format(self.cfg.get('commands', 'skewer_wrapper_cmd'), qc_table_path, 0.75, fq1, fq2, skewer_path, self.cfg.get('commands', 'skewer_cmd'))

        fq0 = ""
        if run_qc:
            fq0 = os.path.join(skewer_path, "out_file-trimmed.fastq.gz")
            fq1 = os.path.join(skewer_path, "out_file-trimmed-pair1.fastq.gz")
            fq2 = os.path.join(skewer_path, "out_file-trimmed-pair2.fastq.gz")
            tools.insert(0, "QC")

        # (0) Readfilter
        cmd_star_filter = "{0} --genomeDir {1} --outFileNamePrefix {2} --readFilesCommand zcat --readFilesIn {3} {4} --outFilterMultimapNmax 100 --outSAMmultNmax 1 --chimSegmentMin 10 --chimJunctionOverhangMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --seedSearchStartLmax 20 --winAnchorMultimapNmax 50 --outSAMtype BAM Unsorted --chimOutType Junctions WithinBAM --outSAMunmapped Within KeepPairs --runThreadN waiting_for_cpu_number".format(self.cfg.get('commands', 'star_cmd'), star_index_path, os.path.join(filtered_reads_path, "{}_".format(sample_id)), fq1, fq2)
        cmd_read_filter = "{0} --input {1}_Aligned.out.bam --output {1}_Aligned.out.filtered.bam".format(self.cfg.get('commands', 'readfilter_cmd'), os.path.join(filtered_reads_path, sample_id))
        # re-define fastq's if filtering is on
#        fq0 = ""
        if filter_reads:
            fq0 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(".fastq.gz", "_filtered_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered.fastq.gz"))
            if run_qc:
                tools.insert(1, "Readfilter")
            else:
                tools.insert(0, "Readfilter")
        cmd_bam_to_fastq = "{0} fastq -0 {1} -1 {2} -2 {3} --threads waiting_for_cpu_number {4}_Aligned.out.filtered.bam".format(self.cfg.get('commands', 'samtools_cmd'), fq0, fq1, fq2, os.path.join(filtered_reads_path, sample_id))
        # (1) Kallisto expression quantification (required for pizzly)
        cmd_kallisto = "{0} quant --threads waiting_for_cpu_number --genomebam --gtf {1} --chromosomes {2} --index {3} --fusion --output-dir waiting_for_output_string {4} {5}".format(self.cfg.get('commands', 'kallisto_cmd'), genes_gtf_path, genome_sizes_path, kallisto_index_path, fq1, fq2)
        # (2) Star expression quantification (required for starfusion and starchip)
        cmd_star = "{0} --genomeDir {1} --outFileNamePrefix waiting_for_output_string --runThreadN waiting_for_cpu_number --runMode alignReads --readFilesIn {2} {3} --readFilesCommand zcat --chimSegmentMin 10 --chimJunctionOverhangMin 10 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --seedSearchStartLmax 20 --winAnchorMultimapNmax 50 --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold".format(self.cfg.get('commands', 'star_cmd'), star_index_path, fq1, fq2)
        # (3) Mapslice
        cmd_extr_fastq1 = "gunzip -c {0} > {1}".format(fq1, fq1[:-3])
        cmd_extr_fastq2 = "gunzip -c {0} > {1}".format(fq2, fq2[:-3])
        cmd_mapsplice = "{0} --chromosome-dir {1} -x {2} -1 {3} -2 {4} --threads waiting_for_cpu_number --output {5} --qual-scale phred33 --bam --seglen 20 --min-map-len 40 --gene-gtf {6} --fusion".format(self.cfg.get('commands', 'mapsplice_cmd'), genome_chrs_path, bowtie_index_path, fq1[:-3], fq2[:-3], mapsplice_path, genes_gtf_path)
        # (4) Fusiocatcher
        cmd_fusioncatcher = "{0} --input {1} --output {2} -p waiting_for_cpu_number".format(self.cfg.get('commands', 'fusioncatcher_cmd'), ",".join([fq1, fq2]), fusioncatcher_path)
        # star-fusion and star-chip can be run upon a previous star run (this is MUST NOT be the star_filter run, but the star_expression run)
        # (5)
        cmd_starfusion = "{0} --chimeric_junction {1} --genome_lib_dir {2} --CPU waiting_for_cpu_number --output_dir {3}".format(self.cfg.get('commands', 'starfusion_cmd'), "{}_Chimeric.out.junction".format(os.path.join(star_path, sample_id)), starfusion_index_path, starfusion_path)
        # (7)
        cmd_starchip = "{0} {1} {2} {3}".format(self.cfg.get("commands", "starchip_cmd"), os.path.join(starchip_path, "starchip"), "{}_Chimeric.out.junction".format(os.path.join(star_path, sample_id)), starchip_param_path)
        # (6) Pizzly
        cmd_pizzly = "{0} -k 29 --gtf {1} --cache {2} --fasta {3} --output {4} {5}".format(self.cfg.get('commands', 'pizzly_cmd'), genes_gtf_path, pizzly_cache_path, genes_fasta_path, os.path.join(pizzly_path, "kallizzy"), os.path.join(kallisto_path, "fusion.txt"))
        cmd_pizzly2 = "{0} {1} {2}".format(self.cfg.get('commands', 'pizzly_cmd2'), "{}.json".format(os.path.join(pizzly_path, "kallizzy")), "{}.json.txt".format(os.path.join(pizzly_path, "kallizzy")))
        # (8) Infusion
        cmd_infusion = "{0} -1 {1} -2 {2} --skip-finished --min-unique-alignment-rate 0 --min-unique-split-reads 0 --allow-non-coding --out-dir {3} {4}".format(self.cfg.get("commands", "infusion_cmd"), fq1, fq2, infusion_path, infusion_cfg_path)
        # (x) Soapfuse
        cmd_soapfuse = "{} -c {} -q {} -i {} -o {} -g {}".format(self.cfg.get('commands', 'soapfuse_wrapper_cmd'), self.cfg.get_path(), qc_table_path, " ".join([fq1, fq2]), soapfuse_path, ref_genome)
        # (9) Data collection
        cmd_fetchdata = "{0} -i {1} -o {2} -s {3} -c {4} --fq1 {5} --fq2 {6} --fusion_support {7} --icam_run".format(self.cfg.get('commands', 'fetchdata_cmd'), output_results_path, fetchdata_path, sample_id, self.cfg.get_path(), fq1, fq2, tool_num_cutoff)
        # (10) De novo assembly of fusion transcripts
        # urla: This is currently still under active development and has not been tested thoroughly
        cmd_denovoassembly = "{0} -i waiting_for_gene_list_input -c {1} -b {2}_Aligned.out.bam -g {3} -t {4} -o waiting_for_assembly_out_dir".format(self.cfg.get('commands', 'denovoassembly_cmd'), self.cfg.get_path(), os.path.join(filtered_reads_path, sample_id), ref_genome, ref_trans)
        # (X) Sample monitoring
        cmd_samples = "{0} --output-file={1} --sample-id={2} --action=append_state --state='".format(self.cfg.get('commands', 'samples_cmd'), self.samples.infile, sample_id)

        # set final lists of executable tools and path
        exe_tools = [
            "QC", #0
            "Readfilter", #1
            "Kallisto", #2
            "Star", #3
            "Mapsplice", #4
            "Fusioncatcher", #5
            "Starfusion", #6
            "Pizzly", #7
            "Starchip", #8
            "Infusion", #9
            "Soapfuse", #10
            "Fetchdata", #11
            "Assembly" #12
            ]
        exe_cmds = [
            " && ".join([cmd_fastqc, cmd_qc_parser, cmd_skewer]), #0
            " && ".join([cmd_star_filter, cmd_read_filter, cmd_bam_to_fastq]), #1
            cmd_kallisto, #2
            cmd_star, #3
            " && ".join([cmd_extr_fastq1, cmd_extr_fastq2, cmd_mapsplice]), #4
            cmd_fusioncatcher, #5
            cmd_starfusion, #6
            " && ".join([cmd_pizzly, cmd_pizzly2]), #7
            cmd_starchip, #8
            cmd_infusion, #9
            cmd_soapfuse, #10
            cmd_fetchdata, #11
            cmd_denovoassembly #12
            ]
        exe_path = [
            qc_path, #0
            filtered_reads_path, #1
            kallisto_path, #2
            star_path, #3
            mapsplice_path, #4
            fusioncatcher_path, #5
            starfusion_path, #6
            pizzly_path, #7
            starchip_path, #8
            infusion_path, #9
            soapfuse_path, #10
            fetchdata_path, #11
            "" #12
            ]

        # create and submit slurm job if the tool is requested and hasn't been run before
        for i, tool in enumerate(exe_tools, 0):
            if not icam_run:
                exe_cmds[i] = exe_cmds[i].replace(" --icam_run", "")
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
                    dependency = Queueing.get_jobs_by_name("Kallisto-{}".format(sample_id))
                elif tool == "Starfusion" or tool == "Starchip":
                    dependency = Queueing.get_jobs_by_name("Star-{}".format(sample_id))
                elif tool == "Fetchdata":
                    dependency = Queueing.get_jobs_by_name(sample_id)
                elif tool == "Assembly":
                    dependency = Queueing.get_jobs_by_name("Fetchdata-{}".format(sample_id))
                elif tool == "ReadFilter":
                    dependency = Queueing.get_jobs_by_name("QC-{}".format(sample_id))
                if filter_reads:
                    dependency.extend(Queueing.get_jobs_by_name("Readfilter-{}".format(sample_id)))
                if run_qc:
                    dependency.extend(Queueing.get_jobs_by_name("QC-{}".format(sample_id)))
                self.logger.debug("Submitting slurm job: CMD - {0}; PATH - {1}; DEPS - {2}".format(cmd, exe_path[i], dependency))
                self.submit_job(uid, cmd, cpu, mem, exe_path[i], dependency)
            else:
                self.logger.info("Skipping {0} as it is not selected for execution (Selected are: {1})".format(tool, tools))

    def submit_job(self, uid, cmd, cores, mem_usage, output_results_folder, dependencies):
        """Submit job to slurm scheduling"""
        # pylint: disable=too-many-arguments
        #         all arguments are required for proper resource management
        already_running = Queueing.get_jobs_by_name("-".join(uid.split("-")[0:2]))
        if not already_running:
            Queueing.submit(uid, cmd, cores, mem_usage, output_results_folder, dependencies, self.partitions, self.userid)
            time.sleep(3)
        else:
            self.logger.error("A job with this application/sample combination is currently running. Skipping {} in order to avoid unintended data loss.".format(uid))

def main():
    """Parse command line arguments and start script"""
    parser = argparse.ArgumentParser(description='Processing of demultiplexed FASTQs')
    # required arguments
    parser.add_argument('-i', '--input', dest='input', help='Specify the fastq folder(s) or fastq file(s) to process.', required=True)
    parser.add_argument('-o', '--output-folder', dest='output_folder', help='Specify the folder to save the results into.', required=True)
    parser.add_argument('-u', '--userid', dest='userid', help='User identifier used for slurm and timing', required=True)
    parser.add_argument('-c', '--config', dest='config', help='Specify config file.', required=True)
    # optional arguments
    parser.add_argument('-p', '--partitions', dest='partitions', help='Comma-separated list of partitions to use for queueing.', default='allNodes')
    parser.add_argument('--tool_support', dest='tool_support', help='The number of tools which are required to support a fusion event.', default="1")
    parser.add_argument('--no_read_filter', dest='filter_reads', help='Run input read data filtering.', default=True, action='store_false')
    parser.add_argument('--no_qc', dest='run_qc', help='Run input read data qc.', default=True, action='store_false')
    parser.add_argument('--version', dest='version', help='Get version info')
    # hidden (not visible in the help display) arguments
    parser.add_argument('--icam_run', dest='icam_run', help=argparse.SUPPRESS, default=False, action='store_true')
    parser.add_argument('--icam_tree', dest='icam_tree', help=argparse.SUPPRESS, default=False, action='store_true')
    parser.add_argument('--overwrite', dest='overwrite', help=argparse.SUPPRESS, default=False, action='store_true')
    args = parser.parse_args()

    # if version is request, print it and exit
    if args.version:
        print(eav.__version__)
        sys.exit(0)

    # create a local copy of the config file in the working dir folder in order to be able to re-run the script
    # urla: todo - the file permission should be set to read only after it was copied to the project folder (doesn't work for icam as "bfx" is a windows file system)

    print(args.run_qc, args.filter_reads)
    proc = Processing(args.config, args.input, args.output_folder, args.partitions, args.userid, args.overwrite)
    proc.run(args.tool_support, args.run_qc, args.filter_reads, args.icam_run, args.icam_tree)

if __name__ == '__main__':
    main()
