#!/usr/bin/env python3

"""
Read input bam file generated by star
Star must have been run with "--outSAMmultNmax 1" and "--chimSegmentMin X"
Aligned pairs are filtered according to their likely of being chimeric with 4 subsequent filter:
    1) the pair is defined as chimeric by star
    2) at least one of the reads in the pair is unmapped
    3) the pair is not a proper pair (unexpected distance)
    4) 10 or more bases of at least one of the reads in the pair are soft-clipped
Multi-mappers are currently discarded as long as they are not belonging to filter "1"
All pairs that are not either "1-4", are defined as unlikely chimeric
A log file is generated after successful completion which contains information about:
    1) number of input vs filtered reads
    2) runtime/speed of the processing
    3) filter counts
    4) 2 qc tests for fast checking that counted pairs are as expected

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""

import sys
import time

import logzero
import pysam
from logzero import logger


class FusionReadFilter(object):
    """Select alignments belonging to putative fusions from an s/bam file"""

    def __init__(self, bam, output):
        """Parameter initialization"""
        self.bam_file = bam
        self.in_bam = pysam.AlignmentFile(bam, "rb")
        self.filtered_bam = pysam.AlignmentFile(output, "wb", template=self.in_bam)
        # a list of ints is used to count pairs assigned to different filters
        # 0: count_input_alignments
        # 1: count_input_pairs
        # 2: count_filtered_pairs
        # 3: count_multimapped
        # 4: count_star_chimeric_alignments
        # 5: count_qcd_alignments
        # 6: count_unmapped
        # 7: count_10bp_s_clip
        # 8: count_proper_pair
        self.counter = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.last_time = 0
        logzero.logfile("{}.fusionReadFilterLog".format(output))

    def classify_pair(self, read1, read2, read1_flag, read2_flag, member_count):
        """Classify a read pairs into different groups"""
        self.counter[1] += 1

        if read1 is None or read2 is None:
            return False

        # count star chimeric
        if member_count > 2:
            self.counter[4] += 1
            return True
        # count multimapping reads
        if read1.get_tag("NH") > 1 or read2.get_tag("NH") > 1:
            self.counter[3] += 1
            return False
        # count pairs containing unmapped mates
        if read1_flag & 0x4 or read2_flag & 0x4:
            self.counter[6] += 1
            return True
        # check if the pair is a "proper pair" (both mapped, adequat distance)
        if not read1_flag & 0x2 or not read2_flag & 0x2:
            self.counter[8] += 1
            return True
        # check if 10 (default) or more soft-clippings (S) are in the alignment
        if (
            read1.query_alignment_length + 10 <= read1.query_length
            or read2.query_alignment_length + 10 <= read2.query_length
        ):
            self.counter[7] += 1
            return True
        return False

    def run(self):
        """Walk linewise through a s/bam file and send proper read pairs to classification"""
        logger.info("Starting fusion read filtering")
        read1 = read2 = None
        read1_flag = read2_flag = None
        self.last_time = time.time()
        last_query = ""
        count_current_query_member = 0
        # read filtering works as follows:
        # urla: I wrote the processing in a way, that implementing the analyses of multi-mapping reads should be straight forward
        #       For this, however, one would need to check in the group of multi-mapping alignments whether there is no possible "normal mapping"
        #       I could so far think of two possible solutions: (1) running "classify_pair" on all possible combinations of the multimapping pairs,
        #       until a "normal mapping" was found. If this is never the case, the pair is written to the filtered out. This, however, will have
        #       a strong effect on runtimes! (2) with the chimeric multi-mapping setting in the most recent star versions. This should actually
        #       work directly as star chimeric classification outranks multi-mapping disposal. It has, however, not been throughly evaluated by the star community

        for read in self.in_bam.fetch(until_eof=True):
            # iterate through alignments as they appear in the file (this is mandatory because
            # (a) we cannot create an index,
            # (b) want to include unmapped reads and
            # (c) have many references in the header during 2nd filtering)
            self.counter[0] += 1
            read_flag = read.flag

            if last_query != read.query_name and self.counter[0] > 1:
                if self.classify_pair(
                    read1, read2, read1_flag, read2_flag, count_current_query_member
                ):
                    self.filtered_bam.write(read1)
                    self.filtered_bam.write(read2)
                    self.counter[2] += 1
                read1 = read2 = None
                read1_flag = read2_flag = None
                count_current_query_member = 0

            count_current_query_member += 1
            # ignore all alignments which are either supplemental, vendor qc'd, secondary or duplicates
            # urla: duplicate flagging/removal prior to running this script would most probably lead to several errors
            #       Nevertheless, it is a very fast check and I would also not recommend that someone does deduplication on rna-seq data!
            if read_flag > 255:
                self.counter[5] += 1
            else:
                # urla: the following commented version should be fine, but the other one is still kept in order to be more error-aware
                #                if read.is_read1:
                #                    read1 = read
                #                else:
                #                    read2 = read
                if not read1 and read.is_read1:
                    read1 = read
                    read1_flag = read_flag
                elif not read2 and read.is_read2:
                    read2 = read
                    read2_flag = read_flag
                else:
                    logger.error(
                        "Neither r1 nor r2??? Read: {0}; R1: {1}; R2: {2}; bamLine: {3}".format(
                            read, read1, read2, self.counter[0]
                        )
                    )
                    sys.exit(1)
            last_query = read.query_name

        # once EOF is reached, the very last pair has to be classified additionally
        if self.classify_pair(
            read1, read2, read1_flag, read2_flag, count_current_query_member
        ):
            self.filtered_bam.write(read1)
            self.filtered_bam.write(read2)
            self.counter[2] += 1

        # close reading/writing stream
        self.in_bam.close()
        self.filtered_bam.close()
        self.print_stats()
        logger.info("Finished fusion read filtering")

    def print_stats(self):
        """print collected statistics to the log file"""
        this_time = time.time()
        time_taken = this_time - self.last_time
        time_taken_1m = 0
        if self.counter[0]:
            time_taken_1m = float(time_taken * 1000000) / self.counter[0]
        self.last_time = this_time

        time_taken_2m = 0
        if self.counter[1]:
            time_taken_2m = float(self.counter[2]) / self.counter[1]

        logger.info(
            "Processed {0} alignments, {1} of {2} pairs remained after filtering ({3:.2%}) ({4:.2f}s / 1M alignments; {5:.2f}s in total)".format(
                self.counter[0],
                self.counter[2],
                self.counter[1],
                time_taken_2m,
                time_taken_1m,
                time_taken,
            )
        )

        qc1 = False
        if self.get_input_read_count_from_star() == self.counter[1]:
            qc1 = True
        qc2 = False
        if (self.counter[4] == self.counter[5]) and (
            (self.counter[0] - self.counter[5]) * 0.5 == self.counter[1]
        ):
            qc2 = True

        # 0: count_input_alignments
        # 1: count_input_pairs
        # 2: count_filtered_pairs
        # 3: count_multimapped
        # 4: count_star_chimeric_alignments
        # 5: count_qcd_alignments
        # 6: count_unmapped
        # 7: count_10bp_s_clip
        # 8: count_proper_pair
        logger.info(
            "Star_chimeric (chim alignment from star):\t{} pairs (filtered)".format(
                self.counter[4]
            )
        )
        logger.info(
            "QC'd (additional Star_chimeric alignment):\t{} alignments (included in above)".format(
                self.counter[5]
            )
        )
        logger.info(
            "Multimapped (1 < x <= 100 equal mappings):\t{} pairs (discarded)".format(
                self.counter[3]
            )
        )
        logger.info(
            "Unmapped (no mapping or >100 multi map):\t{} pairs (filtered)".format(
                self.counter[6]
            )
        )
        logger.info(
            "No proper pair (unexpected read distance):\t{} pairs (filtered)".format(
                self.counter[8]
            )
        )
        logger.info(
            "10bp_s_clip (>9bp soft-clipped in cigar):\t{} pairs (filtered)".format(
                self.counter[7]
            )
        )
        logger.info(
            'Unlikely chimeric ("normal" mappings): \t{} pairs (discarded)'.format(
                self.counter[1]
                - self.counter[4]
                - self.counter[3]
                - self.counter[6]
                - self.counter[8]
                - self.counter[7]
            )
        )
        logger.info("Filter QC1 (fq reads = bam alignments):\t{}".format(qc1))
        logger.info("Filter QC2 (QC'd alignments are chimeric):\t{}".format(qc2))

    def get_input_read_count_from_star(self):
        """Parses a star output log file to get input read counts from the fastq origin"""
        log_file = "{}Log.final.out".format(self.bam_file.rstrip("Aligned.out.bam"))
        with open(log_file, "r") as star_log:
            for line in star_log:
                if line.split("|")[0].strip() == "Number of input reads":
                    return int(line.split("|")[1].strip())
        return -1


def add_read_filter_args(parser):
    parser.add_argument("-i", "--input", dest="input", help="Specify input BAM file")
    parser.add_argument("-o", "--output", dest="output", help="Specify output prefix")
    parser.set_defaults(func=read_filter_command)


def read_filter_command(args):
    read_filter = FusionReadFilter(args.input, args.output)
    read_filter.run()
