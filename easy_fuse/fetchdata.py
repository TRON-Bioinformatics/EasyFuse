#!/usr/bin/env python3

"""
Collect fusion prediction results,
annotate fusions,
run requantification and
combine all results

@author: Tron (PASO), BNT (URLA)
@version: 20190118
"""
import math
import os
import os.path
import subprocess
import sys
from argparse import ArgumentParser

import logzero
from logzero import logger

import easy_fuse.misc.io_methods as IOMethods
from easy_fuse.misc import queueing
from easy_fuse.misc.config import EasyFuseConfiguration
from easy_fuse.misc.samples import SamplesDB


class Fetching(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""

    def __init__(self, sample_path, fetchdata_path, sample_id, config: EasyFuseConfiguration):
        """Parameter initiation and work folder creation."""
        self.sample_path = sample_path
        self.fetchdata_path = fetchdata_path
        self.sample_id = sample_id
        self.samples = SamplesDB(os.path.join(sample_path, os.path.pardir, "samples.db"))
        self.log_file = os.path.join(self.fetchdata_path, "fetchdata.log")
        logzero.logfile(self.log_file)
        self.cfg = config

    def get_pseudo_genome_adjustments_for_star(self, num_len_file):
        """Return the genome size of an associated fasta file calculated by urla_GetFusionSequence_latest.R"""
        with open(num_len_file) as lfile:
            seq_num = int(lfile.readline().rstrip())
            genome_size = int(lfile.readline().rstrip())
        star_genome_chr_bin_n_bits = min(18, int(math.log(genome_size / seq_num, 2)))
        star_genome_sa_index_n_bases = min(14, int(math.log(genome_size, 2) / 2 - 1)) - 2
        logger.debug(
            "Custom genome sequence number: {0} => {1} will be used as bin size parameter for genome storage".format(
                seq_num, star_genome_chr_bin_n_bits))
        logger.debug(
            "Custom Genome Size: {0} bp => {1} will be used as length parameter for SA pre-indexing".format(genome_size,
                                                                                                            star_genome_sa_index_n_bases))
        return (str(star_genome_chr_bin_n_bits), str(star_genome_sa_index_n_bases))

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
        ref_genome = self.cfg["general"]["ref_genome_build"]
        genome_fasta_path = self.cfg["references"]["genome_fasta"]
        genes_adb_path = self.cfg["references"]["genes_adb"]
        genes_tsl_path = self.cfg["references"]["genes_tsl"]

        detected_fusions_path = os.path.join(self.fetchdata_path, "fetched_fusions")
        detected_fusions_file = os.path.join(detected_fusions_path, "Detected_Fusions.csv")
        context_seq_path = os.path.join(self.fetchdata_path, "fetched_contextseqs")
        context_seq_file = os.path.join(context_seq_path, "Context_Seqs.csv")
        filtered_reads_path = os.path.join(self.sample_path, "filtered_reads")
        star_genome_path = os.path.join(context_seq_path, "STAR_idx")
        star_align_path = os.path.join(context_seq_path, "STAR_align")
        star_align_file = os.path.join(star_align_path, "{}_".format(self.sample_id))
        classification_path = os.path.join(self.fetchdata_path, "classification")
        classification_file = os.path.join(classification_path, "classification")

        for folder in [
                self.fetchdata_path,
                detected_fusions_path, 
                context_seq_path,
                star_genome_path, 
                star_align_path,
                classification_path]:
            IOMethods.create_folder(folder)

        # processing steps to perform
        tools = self.cfg["general"]["fd_tools"].split(",")
        fusion_tools = self.cfg["general"]["tools"].split(",")
        cmds = self.cfg["commands"]
        # In case of a liftover, some reference and path must be changed accordingly
        cmd_contextseq_org = ""
        if "liftover" in tools:
            tools.insert(2, "contextSeqBak")
            # for read grepping, we need the original reference on which the first mapping was performed
            cmd_contextseq_org = "fusionannotation " \
                                 "--detected_fusions {0}.bak " \
                                 "--annotation_db {1} " \
                                 "--out_csv {2}.bak " \
                                 "--genome_fasta {3} " \
                                 "--tsl_info {4} " \
                                 "--cis_near_dist {5} " \
                                 "--context_seq_len {6} " \
                                 "--tsl_filter_level {7}".format(
                detected_fusions_file,
                genes_adb_path,
                context_seq_file,
                genome_fasta_path,
                genes_tsl_path,
                self.cfg["general"]["cis_near_distance"],
                self.cfg["general"]["context_seq_len"],
                self.cfg["general"]["tsl_filter"])

            # now, references need to be updated according to the target liftover
            crossmap_chain = self.cfg["liftover"]["crossmap_chain"]
            ref_genome_dest = os.path.basename(crossmap_chain).replace(".", "_").split("_")[2].lower()
            logger.debug(
                "Creating a copy of the detected fusions file due to selection of liftover. Old ({0}) data will be "
                "kept in \"{1}.bak\"".format(ref_genome, detected_fusions_file))
            genome_fasta_path = self.cfg["references"]["genome_fasta_hg37"]
            genes_adb_path = self.cfg["references"]["genes_adb_hg37"]

        # urla - note: tmp hack to get original star input reads for normalization
        with open(os.path.join(classification_path, "Star_org_input_reads.txt"), "w") as infile:
            read_count = self.get_input_read_count_from_star(
                os.path.join(filtered_reads_path, "{}_Aligned.out.bam".format(self.sample_id)))
            if read_count == -1:
                read_count = self.get_input_read_count_from_fastq(fq1)
            infile.write(str(read_count))

        # Define cmd strings for each program
        cmd_fusiondata = "fusiontoolparser " \
                         "-i {0} " \
                         "-o {1} " \
                         "-s {2} " \
                         "-t {3} " \
                         "-f {4} " \
                         "-l {5}".format(
            self.sample_path,
            detected_fusions_path,
            self.sample_id,
            fusion_support,
            self.cfg["general"]["fusiontools"],
            self.log_file)

        cmd_liftover = "liftover " \
                       "-i {0} " \
                       "-l {1} " \
                       "-c {2}".format(
            detected_fusions_file,
            self.log_file,
            self.cfg.config_file)

        cmd_contextseq = "fusionannotation " \
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
            genome_fasta_path,
            genes_tsl_path,
            self.cfg["general"]["cis_near_distance"],
            self.cfg["general"]["context_seq_len"],
            self.cfg["general"]["tsl_filter"])

        cpu = self.cfg["resources"]["fetchdata"].split(",")[0]
        #        mem = cfg.resources["fetchdata"]["mem"]

        cmd_starindex = "{0} " \
                        "--runMode genomeGenerate " \
                        "--runThreadN {1} " \
                        "--limitGenomeGenerateRAM 48000000000 " \
                        "--genomeChrBinNbits waiting_for_bin_size_input " \
                        "--genomeSAindexNbases waiting_for_sa_idx_input " \
                        "--genomeDir {2} --genomeFastaFiles {3}".format(
            cmds["star"],
            cpu,
            star_genome_path,
            "{0}{1}".format(context_seq_file, ".fasta"))

        cmd_staralign_fltr = "{0} " \
                             "--genomeDir {1} " \
                             "--readFilesCommand zcat " \
                             "--readFilesIn {2} {3} " \
                             "--outSAMtype BAM SortedByCoordinate " \
                             "--outFilterMultimapNmax -1 " \
                             "--outSAMattributes Standard " \
                             "--outSAMunmapped None " \
                             "--outFilterMismatchNoverLmax 0.02 " \
                             "--runThreadN {4} " \
                             "--outFileNamePrefix {5}fltr_ " \
                             "--limitBAMsortRAM 48000000000".format(
            cmds["star"],
            star_genome_path,
            fq1,
            fq2,
            cpu,
            star_align_file)

        cmd_bamindex_fltr = "{0} index {1}fltr_Aligned.sortedByCoord.out.bam".format(
            cmds["samtools"],
            star_align_file)

        cmd_requantify_fltr = "requantify " \
                              "-i {0}fltr_Aligned.sortedByCoord.out.bam " \
                              "-o {1}_fltr.tdt " \
                              "-d 10".format(
            star_align_file,
            classification_file)

        (fq1, fq2) = self.samples.get_fastq_files(self.sample_id)

        cmd_staralign_org = "{0} " \
                            "--genomeDir {1} " \
                            "--readFilesCommand zcat " \
                            "--readFilesIn {2} {3} " \
                            "--outSAMtype BAM SortedByCoordinate " \
                            "--outFilterMultimapNmax -1 " \
                            "--outSAMattributes Standard " \
                            "--outSAMunmapped None " \
                            "--outFilterMismatchNoverLmax 0.02 " \
                            "--runThreadN {4} " \
                            "--outFileNamePrefix {5}org_ " \
                            "--limitBAMsortRAM 48000000000".format(
            cmds["star"],
            star_genome_path,
            fq1,
            fq2,
            cpu,
            star_align_file)

        cmd_bamindex_org = "{0} index {1}org_Aligned.sortedByCoord.out.bam".format(
            cmds["samtools"],
            star_align_file)

        cmd_requantify_org = "requantify " \
                             "-i {0}org_Aligned.sortedByCoord.out.bam " \
                             "-o {1}_org.tdt -d 10".format(
            star_align_file,
            classification_file)

        # for testing, based on debug. should be removed if merged to original
        cmd_read_filter2 = "getRequantReads " \
                           "--input {0}_Aligned.out.bam " \
                           "--input2 {1}.debug " \
                           "--output {0}_Aligned.out.filtered2.bam".format(
            os.path.join(filtered_reads_path, self.sample_id),
            context_seq_file)

        # re-define fastq's if filtering is on (default)
        fq0 = ""
        if "readfilter" in fusion_tools:
            fq0 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace("R1", "R0").replace(".fastq.gz",
                                                                                                      "_filtered2_singles.fastq.gz"))
            fq1 = os.path.join(filtered_reads_path, os.path.basename(fq1).replace(".fastq.gz", "_filtered2.fastq.gz"))
            fq2 = os.path.join(filtered_reads_path, os.path.basename(fq2).replace(".fastq.gz", "_filtered2.fastq.gz"))

        cmd_bam_to_fastq = "{0} fastq " \
                           "-0 {1} " \
                           "-1 {2} " \
                           "-2 {3} " \
                           "--threads {5} {4}_Aligned.out.filtered2.bam".format(
            cmds["samtools"],
            fq0,
            fq1,
            fq2,
            os.path.join(filtered_reads_path, self.sample_id),
            cpu)

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
                             "--runThreadN {4} " \
                             "--outFileNamePrefix {5}best_ " \
                             "--limitBAMsortRAM 48000000000".format(
            cmds["star"],
            star_genome_path,
            fq1,
            fq2,
            cpu,
            star_align_file)

        cmd_bamindex_best = "{0} index {1}best_Aligned.sortedByCoord.out.bam".format(
            cmds["samtools"],
            star_align_file)

        cmd_requantify_best = "requantify " \
                              "-i {0}best_Aligned.sortedByCoord.out.bam " \
                              "-o {1}_best.tdt -d 10".format(
            star_align_file,
            classification_file)

        # set final lists of executable tools and path
        exe_tools = [
            "fusiongrep", #1
            "liftover", #2
            "contextSeqBak",
            "contextseq", #3
            "starindex", #4
            "staralignFltr", #5
            "bamindexFltr", #6
            "requantifyFltr", #7
            "staralignOrg", #8
            "bamindexOrg", #9
            "requantifyOrg", #10
            "readFilter2", #11
            "readFilter2b", #12
            "staralignBest", #13
            "bamindexBest", #14
            "requantifyBest" #15
            ]
        exe_cmds = [
            cmd_fusiondata,  # 1
            cmd_liftover,  # 2
            cmd_contextseq_org,
            cmd_contextseq,  # 3
            cmd_starindex,  # 4
            cmd_staralign_fltr,  # 5
            cmd_bamindex_fltr,  # 6
            cmd_requantify_fltr,  # 7
            cmd_staralign_org,  # 8
            cmd_bamindex_org,  # 9
            cmd_requantify_org,  # 10
            cmd_read_filter2,  # 11
            cmd_bam_to_fastq,  # 12
            cmd_staralign_best,  # 13
            cmd_bamindex_best,  # 14
            cmd_requantify_best  # 15
        ]
        exe_dependencies = [
            "",  # 1
            detected_fusions_file,  # 2
            detected_fusions_file,
            detected_fusions_file,  # 3
            "{0}{1}".format(context_seq_file, ".fasta.info"),  # 4
            star_genome_path,  # 5
            "{}fltr_Aligned.sortedByCoord.out.bam".format(star_align_file),  # 6
            "",  # 7
            star_genome_path,  # 8
            "{}org_Aligned.sortedByCoord.out.bam".format(star_align_file),  # 9
            "",  # 10
            "",  # 11
            "",  # 12
            star_genome_path,  # 13
            "{}best_Aligned.sortedByCoord.out.bam".format(star_align_file),  # 14
            ""  # 15
        ]

        # create and submit slurm job if the tool is requested and hasn't been run before
        module_file = self.cfg["general"]["build_env"]
        for i, tool in enumerate(exe_tools, 0):
            if tool in tools:
                if not exe_dependencies[i] or os.path.exists(exe_dependencies[i]):
                    logger.info("Starting {}".format(tool))
                    if tool == "starindex":  # the genome size required for the genomeSAindexNbases parameter is not known before now
                        (star_bin, star_sa) = self.get_pseudo_genome_adjustments_for_star(
                            "{0}{1}".format(context_seq_file, ".fasta.info"))
                        exe_cmds[i] = exe_cmds[i].replace("waiting_for_bin_size_input", star_bin)
                        exe_cmds[i] = exe_cmds[i].replace("waiting_for_sa_idx_input", star_sa)
                    logger.debug("Executing: {}".format(exe_cmds[i]))
                    queueing.submit("", exe_cmds[i].split(" "), "", "", "", "", "", "", "", "", module_file, "none")
                else:
                    logger.error(
                        "Could not run {0} due to the missing dependency {1}".format(tool, exe_dependencies[i]))
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
    parser.add_argument('-c', '--config-file', dest='config_file', required=True,
                        help='Specify alternative config file to use for your analysis')
    parser.add_argument('--fq1', dest='fq1', help='Input read1 file for re-quantification', default="", required=False)
    parser.add_argument('--fq2', dest='fq2', help='Input read2 file for re-quantification', default="", required=False)
    parser.add_argument('--fusion_support', dest='fusion_support',
                        help='Number of fusion tools which must predict the fusion event', default=1)
    args = parser.parse_args()

    config = EasyFuseConfiguration(args.config_file)
    # processing
    # 1: fusion tool parser
    # 2: get context seqs
    # 3: re-quantify
    # 4: get expression?
    # ...
    for i in range(1, int(args.fusion_support) + 1):
        proc = Fetching(args.input, args.output, args.sample, config)
        proc.run(i, args.fq1, args.fq2)


if __name__ == '__main__':
    main()
