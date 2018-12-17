"""
Collect fusion prediction results,
annotate fusions,
run requantification and
combine all results

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""

from __future__ import print_function
import sys
import os
import os.path
import time
import math
from argparse import ArgumentParser
import misc.queue as Queueing
from misc.config import Config
from misc.samples import Samples
from misc.logger import Logger
import misc.io_methods as IOMethods

# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class Fetching(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""
    def __init__(self, scratch_path, fetchdata_path, sample_id, config):
        """Parameter initiation and work folder creation."""
        self.cfg = config
        self.scratch_path = scratch_path
        self.fetchdata_path = fetchdata_path
        self.sample_id = sample_id
        #self.tools = Samples(os.path.join(scratch_path, os.path.pardir, os.path.pardir, "samples.csv")).get_tool_list_from_state(self.sample_id)
        self.sample = Samples(os.path.join(scratch_path, os.path.pardir, os.path.pardir, "samples.csv"))

        self.logger = Logger(os.path.join(self.fetchdata_path, "fetchdata_" + str(int(round(time.time()))) + ".log"))

    def get_pseudo_genome_adjustments_for_star(self, num_len_file): # wrong pylint error due to long name => pylint: disable=C0103
        """Return the genome size of an associated fasta file calculated by urla_GetFusionSequence_latest.R"""
        seq_num = 0
        genome_size = 0
        with open(num_len_file) as lfile:
            seq_num = int(lfile.next())
            genome_size = int(lfile.next())
        star_genome_chr_bin_n_bits = min(18, int(math.log(genome_size / seq_num, 2)))
        star_genome_sa_index_n_bases = min(14, int(math.log(genome_size, 2) / 2 - 1))
        self.logger.debug("Custom genome sequence number: {0} => {1} will be used as bin size parameter for genome storage".format(seq_num, star_genome_chr_bin_n_bits))
        self.logger.debug("Custom Genome Size: {0} bp => {1} will be used as length parameter for SA pre-indexing".format(genome_size, star_genome_sa_index_n_bases))
        return(str(star_genome_chr_bin_n_bits), str(star_genome_sa_index_n_bases))

    def run(self, ref_genome, ref_trans, fusion_support, fq1, fq2):
        """Identification of fastq files and initiation of processing"""
        # print sample id
        # execute processing pipe
        # sampleID = ...
        self.logger.info("Fetching in sample {}".format(self.sample_id))
        if not fq1 or not fq2:
            self.logger.debug("Either ReadFile 1 or 2 or both are missing, trying to get original files from samples.csv")
            (fq1, fq2) = self.sample.get_fastq_files(self.sample_id)
        self.execute_pipeline((fq1, fq2), ref_genome, ref_trans, fusion_support)

    # urla: there are a lot of local variables declarated in the following method.
    #       Although this could be reduced quite strongly, readability would be strongly reduced as well
    #       pylint:disable=R0914
    def execute_pipeline(self, (fq1, fq2), ref_genome, ref_trans, fusion_support):
        """Create sample specific subfolder structuge and run tools on fastq files"""

		 # Genome/Gene references to use
        genome_fastadir_path = self.cfg.get('references', ref_trans + '_genome_fastadir_' + ref_genome)
        genes_csv_path = self.cfg.get('references', ref_trans + '_genes_csv_' + ref_genome)

        fetchdata_current_path = os.path.join(self.fetchdata_path, "fd_{}_tool".format(fusion_support))
        detected_fusions_path = os.path.join(fetchdata_current_path, "fetched_fusions")
        detected_fusions_file = os.path.join(detected_fusions_path, "Detected_Fusions.csv")
        context_seq_path = os.path.join(fetchdata_current_path, "fetched_contextseqs")
        context_seq_file = os.path.join(context_seq_path, "Context_Seqs.csv")
        star_genome_path = os.path.join(context_seq_path, "STAR_idx")
        star_align_path = os.path.join(context_seq_path, "STAR_align")
        star_align_file = os.path.join(star_align_path, "{}_".format(self.sample_id))
        classification_path = os.path.join(fetchdata_current_path, "classification")
        classification_file = os.path.join(classification_path, "classification.tdt")

        for folder in [
                fetchdata_current_path,
                detected_fusions_path, context_seq_path,
                star_genome_path, star_align_path,
                classification_path
            ]:
            IOMethods.create_folder(folder, self.logger)

        # Define cmd strings for each program
        fusiontools = "Mapsplice,Fusioncatcher,Starfusion,Starchip,Infusion"
        cmd_fusiondata = "{0} -i {1} -o {2} -s {3} -t {4} -f {5} -l {6}".format(self.cfg.get('commands', 'fetch_fusiondata_cmd'), self.scratch_path, detected_fusions_path, self.sample_id, fusion_support, fusiontools, self.logger.get_path())
        cmd_contextseq = "{0} --input_detected_fusions {1} --fasta_genome_dir {2} --ensembl_csv {3} --output {4}".format(self.cfg.get('commands', 'fetch_context_cmd'), detected_fusions_file, genome_fastadir_path, genes_csv_path, context_seq_file)

        cpu = 12
        cmd_starindex = "{0} --runMode genomeGenerate --runThreadN {1} --limitGenomeGenerateRAM 30000000000 --genomeChrBinNbits waiting_for_bin_size_input --genomeSAindexNbases waiting_for_sa_idx_input --genomeDir {2} --genomeFastaFiles {3}".format(self.cfg.get('commands', 'star_cmd'), cpu, star_genome_path, "{0}{1}".format(context_seq_file, ".fasta"))
        cmd_staralign = "{0} --genomeDir {1} --readFilesCommand zcat --readFilesIn {2} {3} --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax -1 --outSAMattributes Standard --outSAMunmapped None --outFilterMismatchNoverLmax 0.02 --runThreadN {4} --outFileNamePrefix {5} --limitBAMsortRAM 30000000000".format(self.cfg.get('commands', 'star_cmd'), star_genome_path, fq1, fq2, cpu, star_align_file)
        cmd_bamindex = "{0} index {1}Aligned.sortedByCoord.out.bam".format(self.cfg.get('commands', 'samtools_cmd'), star_align_file)
        cmd_requantify = "{0} -i {1}Aligned.sortedByCoord.out.bam -o {2} -d 3".format(self.cfg.get('commands', 'classification_cmd'), star_align_file, classification_file)

        tools = self.cfg.get('general', 'fd_tools').split(",")

        # set final lists of executable tools and path
        exe_tools = [
            "Fusiongrep", #1
            "Contextseq", #2
            "Starindex", #3
            "Staralign", #4
            "Bamindex", #5
            "Requantify" #6
            ]
        exe_cmds = [
            cmd_fusiondata, #1
            cmd_contextseq, #2
            cmd_starindex, #3
            cmd_staralign, #4
            cmd_bamindex, #5
            cmd_requantify #6
            ]
        exe_dependencies = [
            "", #1
            detected_fusions_file, #2
            "{0}{1}".format(context_seq_file, ".fasta.info"), #3
            star_genome_path, #4
            "{}Aligned.sortedByCoord.out.bam".format(star_align_file), #5
            "" #6
            ]

        # create and submit slurm job if the tool is requested and hasn't been run before
        for i, tool in enumerate(exe_tools, 0):
            if tool in tools:
                if not exe_dependencies[i] or os.path.exists(exe_dependencies[i]):
                    self.logger.info("Starting {}".format(tool))
                    if tool == "Starindex": # the genome size required for the genomeSAindexNbases parameter is not known before now
                        (star_bin, star_sa) = self.get_pseudo_genome_adjustments_for_star("{0}{1}".format(context_seq_file, ".fasta.info"))
                        exe_cmds[i] = exe_cmds[i].replace("waiting_for_bin_size_input", star_bin)
                        exe_cmds[i] = exe_cmds[i].replace("waiting_for_sa_idx_input", star_sa)
                    #cmd = " && ".join([exe_cmds[i], cmd_samples + uid + "'"])
                    self.logger.debug("Executing: {}".format(exe_cmds[i]))
                    Queueing._submit_nonqueue(exe_cmds[i].split(" "))
                else:
                    self.logger.error("Could not run {0} due to the missing dependency {1}".format(tool, exe_dependencies[i]))
                    sys.exit(1)
            else:
                self.logger.debug("Skipping {0} as it is not selected for execution (Selected are: {1})".format(tool, tools))

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description='Processing of demultiplexed FASTQs')
    # required arguments
    parser.add_argument('-i', '--input', dest='input', help='Data input directory.', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Data output directory.', required=True)
    parser.add_argument('-s', '--sample', dest='sample', help='Specify the sample to process.', required=True)
    parser.add_argument('-t', '--reference-transcripts', dest='ref_trans', choices=['ucsc', 'refseq', 'ensembl'], help='Specify the reference transcript database name.', default='ensembl', required=True)
    parser.add_argument('-g', '--reference-genome', dest='ref_genome', choices=['hg18', 'hg19', 'hg38', 'mm9', 'mm10'], help='Specify the reference genome name.', default='hg38', required=True)
    parser.add_argument('-c', '--config', dest='config', help='Specify config file.', default="", required=True)
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
        proc = Fetching(args.input, args.output, args.sample, Config(args.config))
        proc.run(args.ref_genome, args.ref_trans, i, args.fq1, args.fq2)

if __name__ == '__main__':
    main()