#!/usr/bin/env python3

"""
Collect fusion prediction results,
annotate fusions,
run requantification and
combine all results

@author: Tron (PASO), BNT (URLA)
@version: 20190118
"""

from argparse import ArgumentParser
from configparser import ConfigParser
import json
import math
import os
import os.path
import subprocess
import sys

from easy_fuse.misc import queueing
# custom imports
from easy_fuse.misc.samples import SamplesDB
import easy_fuse.misc.io_methods as IOMethods
import logzero
from logzero import logger


class Fetching(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""
    def __init__(self, scratch_path, fetchdata_path, sample_id, cfg_file):
        """Parameter initiation and work folder creation."""
        self.scratch_path = scratch_path
        self.fetchdata_path = fetchdata_path
        self.sample_id = sample_id
        self.samples = SamplesDB(os.path.join(scratch_path, os.path.pardir, "samples.db"))
        logzero.logfile(os.path.join(self.fetchdata_path, "fetchdata.log"))

        self.cfg = None


        default_cfg_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini")

        if not cfg_file and os.path.exists(default_cfg_path):
            cfg_file = default_cfg_path
        elif not cfg_file and not os.path.exists(default_cfg_path):
            logger.error("Could not find default config file (path={}). Exiting.".format(default_cfg_path))
            sys.exit(1)

        if cfg_file.endswith("ini"):
            self.cfg = ConfigParser()
            self.cfg.read(cfg_file)
        elif cfg_file.endswith("json"):
            with open(cfg_file) as config_file:
                self.cfg = json.load(config_file)


    def get_pseudo_genome_adjustments_for_star(self, num_len_file):
        """Return the genome size of an associated fasta file calculated by urla_GetFusionSequence_latest.R"""
        seq_num = 0
        genome_size = 0
        with open(num_len_file) as lfile:
            seq_num = int(lfile.readline().rstrip())
            genome_size = int(lfile.readline().rstrip())
        star_genome_chr_bin_n_bits = min(18, int(math.log(genome_size / seq_num, 2)))
        star_genome_sa_index_n_bases = min(14, int(math.log(genome_size, 2) / 2 - 1)) - 2
        logger.debug("Custom genome sequence number: {0} => {1} will be used as bin size parameter for genome storage".format(seq_num, star_genome_chr_bin_n_bits))
        logger.debug("Custom Genome Size: {0} bp => {1} will be used as length parameter for SA pre-indexing".format(genome_size, star_genome_sa_index_n_bases))
        return(str(star_genome_chr_bin_n_bits), str(star_genome_sa_index_n_bases))

    @staticmethod
    def get_input_read_count_from_star(star_out_bam):
        """Parses a star output log file to get input read counts from the fastq origin"""
        log_file = "{}Log.final.out".format(star_out_bam.rstrip("Aligned.out.bam"))
        if not os.path.exists(log_file):
            return -1
        with open(log_file, "r") as star_log:
            for line in star_log:
                if line.split("|")[0].strip() == "Number of input reads":
                    return int(line.split("|")[1].strip())
        return -1

    @staticmethod
    def get_input_read_count_from_fastq(fastq):
        """Parses input FASTQ to get read count"""
        ps = subprocess.Popen(("zcat", fastq), stdout=subprocess.PIPE)
        result = subprocess.check_output(("wc", "-l"), stdin=ps.stdout)
        return int(result) / 4

    def run(self, fusion_support, fq1, fq2):
        """Identification of fastq files and initiation of processing"""
        # print sample id
        # execute processing pipe
        # sampleID = ...
        logger.info("Fetching in sample {}".format(self.sample_id))
        if not fq1 or not fq2:
            logger.debug("Either ReadFile 1 or 2 or both are missing, trying to get original files from samples.csv")
            logger.debug(self.sample_id)
            logger.debug(self.samples.db_path)
            (fq1, fq2) = self.samples.get_fastq_files(self.sample_id)
        self.execute_pipeline(fq1, fq2, fusion_support)

    # urla: there are a lot of local variables declarated in the following method.
    #       Although this could be reduced quite strongly, readability would be strongly reduced as well
    #       pylint:disable=R0914
    def execute_pipeline(self, fq1, fq2, fusion_support):
        """Create sample specific subfolder structuge and run tools on fastq files"""

        # Genome/Gene references to use
        ref_trans = self.cfg["general"]["ref_trans_version"]
        ref_genome = self.cfg["general"]["ref_genome_build"]
        genome_fasta_path = self.cfg["references"]["genome_fasta"]
        genes_adb_path = self.cfg["references"]["genes_adb"]
        genes_tsl_path = self.cfg["references"]["genes_tsl"]

        fetchdata_current_path = os.path.join(self.fetchdata_path, "fd_{}_tool".format(fusion_support))
        detected_fusions_path = os.path.join(fetchdata_current_path, "fetched_fusions")
        detected_fusions_file = os.path.join(detected_fusions_path, "Detected_Fusions.csv")
        context_seq_path = os.path.join(fetchdata_current_path, "fetched_contextseqs")
        context_seq_file = os.path.join(context_seq_path, "Context_Seqs.csv")
        filtered_reads_path = os.path.join(self.scratch_path, "filtered_reads")
        star_genome_path = os.path.join(context_seq_path, "STAR_idx")
        star_align_path = os.path.join(context_seq_path, "STAR_align")
        star_align_file = os.path.join(star_align_path, "{}_".format(self.sample_id))
        classification_path = os.path.join(fetchdata_current_path, "classification")
        classification_file = os.path.join(classification_path, "classification")

        for folder in [
                fetchdata_current_path,
                detected_fusions_path, 
                context_seq_path,
                star_genome_path, 
                star_align_path,
                classification_path
            ]:
            IOMethods.create_folder(folder)

        # processing steps to perform
        tools = self.cfg["general"]["fd_tools"].split(",")
        fusion_tools = self.cfg["general"]["tools"].split(",")
        module_dir = self.cfg["general"]["module_dir"]
        cmds = self.cfg["commands"]
        # In case of a liftover, some reference and path must be changed accordingly
        cmd_contextseq_org = ""
        if "Liftover" in tools:
            tools.insert(2, "ContextSeqBak")
            # for read grepping, we need the original reference on which the first mapping was performed
            cmd_contextseq_org = "python {0} --detected_fusions {1}.bak --annotation_db {2} --out_csv {3}.bak --genome_fasta {4} --tsl_info {5} --cis_near_dist {6} --context_seq_len {7} --tsl_filter_level {8}".format(os.path.join(module_dir,
                                                                                                                                                                                                                                      "fusionannotation.py"), detected_fusions_file, genes_adb_path, context_seq_file, genome_fasta_path, genes_tsl_path, self.cfg["general"]["cis_near_distance"], self.cfg["general"]["context_seq_len"], self.cfg["general"]["tsl_filter"])
            # now, references need to be updated according to the target liftover
            crossmap_chain = self.cfg["liftover"]["crossmap_chain"]
            ref_genome_dest = os.path.basename(crossmap_chain).replace(".", "_").split("_")[2].lower()
            logger.debug("Creating a copy of the detected fusions file due to selection of liftover. Old ({0}) data will be kept in \"{1}.bak\"".format(ref_genome, detected_fusions_file))
            genome_fasta_path = self.cfg["references"]["genome_fasta_hg37"]
            genes_adb_path = self.cfg["references"]["genes_adb_hg37"]

        # urla - note: tmp hack to get original star input reads for normalization
        with open(os.path.join(classification_path, "Star_org_input_reads.txt"), "w") as infile:
            read_count = self.get_input_read_count_from_star(os.path.join(filtered_reads_path, "{}_Aligned.out.bam".format(self.sample_id)))
            if read_count == -1:
                read_count = self.get_input_read_count_from_fastq(fq1)
            infile.write(str(read_count))
        # Define cmd strings for each program
        cmd_fusiondata = "{0} -i {1} -o {2} -s {3} -t {4} -f {5} -l {6}".format(os.path.join(module_dir,
                                                                                                 "fusiontoolparser.py"), self.scratch_path, detected_fusions_path, self.sample_id, fusion_support, self.cfg["general"]["fusiontools"], logger.get_path())
        cmd_liftover = "{0} -i {1} -l {2}".format(os.path.join(module_dir, "misc", "liftover.py"), detected_fusions_file, logger.get_path())
        cmd_contextseq = "{0} --detected_fusions {1} --annotation_db {2} --out_csv {3} --genome_fasta {4} --tsl_info {5} --cis_near_dist {6} --context_seq_len {7} --tsl_filter_level {8}".format(os.path.join(module_dir,
                                                                                                                                                                                                               "fusionannotation.py"), detected_fusions_file, genes_adb_path, context_seq_file, genome_fasta_path, genes_tsl_path, self.cfg["general"]["cis_near_distance"], self.cfg["general"]["context_seq_len"], self.cfg["general"]["tsl_filter"])
        cpu = self.cfg["resources"]["fetchdata"].split(",")[0]
#        mem = cfg.resources["fetchdata"]["mem"]
        cmd_starindex = "{0} --runMode genomeGenerate --runThreadN {1} --limitGenomeGenerateRAM 48000000000 --genomeChrBinNbits waiting_for_bin_size_input --genomeSAindexNbases waiting_for_sa_idx_input --genomeDir {2} --genomeFastaFiles {3}".format(cmds["star"], cpu, star_genome_path, "{0}{1}".format(context_seq_file, ".fasta"))
        cmd_staralign_fltr = "{0} --genomeDir {1} --readFilesCommand zcat --readFilesIn {2} {3} --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax -1 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 --runThreadN {4} --outFileNamePrefix {5}fltr_ --limitBAMsortRAM 48000000000".format(cmds["star"], star_genome_path, fq1, fq2, cpu, star_align_file)
        cmd_bamindex_fltr = "{0} index {1}fltr_Aligned.sortedByCoord.out.bam".format(cmds["samtools"], star_align_file)
        cmd_requantify_fltr = "{0} -i {1}fltr_Aligned.sortedByCoord.out.bam -o {2}_fltr.tdt -d 10".format(os.path.join(module_dir,
                                                                                                                       "requantify.py"), star_align_file, classification_file)
        (fq1, fq2) = self.samples.get_fastq_files(self.sample_id)
        cmd_staralign_org = "{0} --genomeDir {1} --readFilesCommand zcat --readFilesIn {2} {3} --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax -1 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 --runThreadN {4} --outFileNamePrefix {5}org_ --limitBAMsortRAM 48000000000".format(cmds["star"], star_genome_path, fq1, fq2, cpu, star_align_file)
        cmd_bamindex_org = "{0} index {1}org_Aligned.sortedByCoord.out.bam".format(cmds["samtools"], star_align_file)
        cmd_requantify_org = "{0} -i {1}org_Aligned.sortedByCoord.out.bam -o {2}_org.tdt -d 10".format(os.path.join(module_dir,
                                                                                                                    "requantify.py"), star_align_file, classification_file)
        # for testing, based on debug. should be removed if merged to original
        cmd_read_filter2 = "{0} --input {1}_Aligned.out.bam --input2 {2}.debug --output {1}_Aligned.out.filtered2.bam".format(os.path.join(module_dir,
                                                                                                                                           "getRequantReads.py"), os.path.join(filtered_reads_path, self.sample_id), context_seq_file)
        # re-define fastq's if filtering is on (default)
        fq0 = ""
        if "Readfilter" in fusion_tools:
            fq0 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(".fastq.gz", "_filtered2_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered2.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered2.fastq.gz"))
        cmd_bam_to_fastq = "{0} fastq -0 {1} -1 {2} -2 {3} --threads {5} {4}_Aligned.out.filtered2.bam".format(cmds["samtools"], fq0, fq1, fq2, os.path.join(filtered_reads_path, self.sample_id), cpu)
        # allow soft-clipping? Specificity? --alignEndsType EndToEnd
        cmd_staralign_best = "{0} --genomeDir {1} --readFilesCommand zcat --readFilesIn {2} {3} --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax -1 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 --runThreadN {4} --outFileNamePrefix {5}best_ --limitBAMsortRAM 48000000000".format(cmds["star"], star_genome_path, fq1, fq2, cpu, star_align_file)
        cmd_bamindex_best = "{0} index {1}best_Aligned.sortedByCoord.out.bam".format(cmds["samtools"], star_align_file)
        cmd_requantify_best = "{0} -i {1}best_Aligned.sortedByCoord.out.bam -o {2}_best.tdt -d 10".format(os.path.join(module_dir,
                                                                                                                       "requantify.py"), star_align_file, classification_file)


        # set final lists of executable tools and path
        exe_tools = [
            "Fusiongrep", #1
            "Liftover", #2
            "ContextSeqBak",
            "Contextseq", #3
            "Starindex", #4
            "StaralignFltr", #5
            "BamindexFltr", #6
            "RequantifyFltr", #7
            "StaralignOrg", #8
            "BamindexOrg", #9
            "RequantifyOrg", #10
            "ReadFilter2", #11
            "ReadFilter2b", #12
            "StaralignBest", #13
            "BamindexBest", #14
            "RequantifyBest" #15
            ]
        exe_cmds = [
            cmd_fusiondata, #1
            cmd_liftover, #2
            cmd_contextseq_org,
            cmd_contextseq, #3
            cmd_starindex, #4
            cmd_staralign_fltr, #5
            cmd_bamindex_fltr, #6
            cmd_requantify_fltr, #7
            cmd_staralign_org, #8
            cmd_bamindex_org, #9
            cmd_requantify_org, #10
            cmd_read_filter2, #11
            cmd_bam_to_fastq, #12
            cmd_staralign_best, #13
            cmd_bamindex_best, #14
            cmd_requantify_best #15
            ]
        exe_dependencies = [
            "", #1
            detected_fusions_file, #2
            detected_fusions_file,
            detected_fusions_file, #3
            "{0}{1}".format(context_seq_file, ".fasta.info"), #4
            star_genome_path, #5
            "{}fltr_Aligned.sortedByCoord.out.bam".format(star_align_file), #6
            "", #7
            star_genome_path, #8
            "{}org_Aligned.sortedByCoord.out.bam".format(star_align_file), #9
            "", #10
            "", #11
            "", #12
            star_genome_path, #13
            "{}best_Aligned.sortedByCoord.out.bam".format(star_align_file), #14
            "" #15
            ]

        # create and submit slurm job if the tool is requested and hasn't been run before
        module_file = os.path.join(module_dir, "build_env.sh")
        for i, tool in enumerate(exe_tools, 0):
            if tool in tools:
                if not exe_dependencies[i] or os.path.exists(exe_dependencies[i]):
                    logger.info("Starting {}".format(tool))
                    if tool == "Starindex": # the genome size required for the genomeSAindexNbases parameter is not known before now
                        (star_bin, star_sa) = self.get_pseudo_genome_adjustments_for_star("{0}{1}".format(context_seq_file, ".fasta.info"))
                        exe_cmds[i] = exe_cmds[i].replace("waiting_for_bin_size_input", star_bin)
                        exe_cmds[i] = exe_cmds[i].replace("waiting_for_sa_idx_input", star_sa)
                    logger.debug("Executing: {}".format(exe_cmds[i]))
                    queueing.submit("", exe_cmds[i].split(" "), "", "", "", "", "", "", "", "", module_file, "none")
                else:
                    logger.error("Could not run {0} due to the missing dependency {1}".format(tool, exe_dependencies[i]))
                    sys.exit(1)
            else:
                logger.debug("Skipping {0} as it is not selected for execution (Selected are: {1})".format(tool, tools))

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description='Processing of demultiplexed FASTQs')
    # required arguments
    parser.add_argument('-i', '--input', dest='input', help='Data input directory.', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Data output directory.', required=True)
    parser.add_argument('-s', '--sample', dest='sample', help='Specify the sample to process.', required=True)
    parser.add_argument('-c', '--config-file', dest='config_file', help='Specify alternative config file to use for your analysis')
    parser.add_argument('--fq1', dest='fq1', help='Input read1 file for re-quantification', default="", required=False)
    parser.add_argument('--fq2', dest='fq2', help='Input read2 file for re-quantification', default="", required=False)
    parser.add_argument('--fusion_support', dest='fusion_support', help='Number of fusion tools which must predict the fusion event', default=1)
    args = parser.parse_args()

    # processing
    # 1: fusion tool parser
    # 2: get context seqs
    # 3: re-quantify
    # 4: get expression?
    # ...
    for i in range(1, int(args.fusion_support) + 1):
        proc = Fetching(args.input, args.output, args.sample, os.path.abspath(args.config_file))
        proc.run(i, args.fq1, args.fq2)

if __name__ == '__main__':
    main()
