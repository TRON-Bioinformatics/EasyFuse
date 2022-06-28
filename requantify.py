#!/usr/bin/env python

"""
Junction/Spanning read pair counting

@author: BNT (URLA)
@version: 20190118
"""

from __future__ import print_function, division
import os.path
from argparse import ArgumentParser
#from urla_logger_latest import Logger
import pysam # pysam is not available for windows (where I run pylint) => pylint: disable=E0401


# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class Requantification(object):
    """Select alignments belonging to putative fusions from an s/bam file"""
    def __init__(self, bam, output, bp_distance_threshold):
        """Parameter initialization"""
        self.in_bam = pysam.AlignmentFile(bam, "rb")
        self.output = output
        self.bp_distance_threshold = int(bp_distance_threshold)
        self.fusion_seq_dict = {}
        #self.logger = Logger("{}.fusionReadFilterLog".format(output))
        self.input_read_count = 0

    def quantify_read_groups(self, read_buffer, reference_base):
        """Count junction/spanning pairs on ft and wt1/2"""
        anchor_ft = 0
        anchor_wt1 = 0
        anchor_wt2 = 0
        # the read_buffer contains all reads mapping to the fusion transcript and the wildtyp backgrounds
        # the read_group is list of reads (incl multimappers) belonging to a single read pair
        # reads from this read group can map independently to ft, wt1 and wt2
        for read_group in read_buffer:
            junction_ft = 0 # 0 = no junction overlapping read, 1+ = at least one of the paired reads is a junction read
            spanning_ft = [0, 0] # left: paired read cnt fully left of the fusion bp, right: paired read cnt fulls right of the fusion bp
            junction_wt1 = 0 # like is_junction_ft but for wt1
            spanning_wt1 = [0, 0] # like is_spanning_ft but for wt1
            junction_wt2 = 0 # like is_junction_ft but for wt2
            spanning_wt2 = [0, 0] # like is_spanning_ft but for wt2

            # each processed read belongs to one of the three defined read_groups
            for read in read_buffer[read_group]:
                if read.reference_name.endswith("ft"):
                    junction_ft, spanning_ft, anchor_ft = self.count_junc_span(read, self.fusion_seq_dict[reference_base][0], junction_ft, spanning_ft, anchor_ft)
                elif read.reference_name.endswith("wt1"):
                    junction_wt1, spanning_wt1, anchor_wt1 = self.count_junc_span(read, self.fusion_seq_dict[reference_base][6], junction_wt1, spanning_wt1, anchor_wt1)
                elif read.reference_name.endswith("wt2"):
                    junction_wt2, spanning_wt2, anchor_wt2 = self.count_junc_span(read, self.fusion_seq_dict[reference_base][12], junction_wt2, spanning_wt2, anchor_wt2)

            #print("stored before: {}".format(self.fusion_seq_dict[reference_base]))
            #print("jft: {}, sft: {}, jwt: {}, swt: {}".format(is_junction_ft, is_spanning_ft, is_junction_wt, is_spanning_wt))
            # update junction and spanning counts
            self.update_counts(reference_base, junction_ft, spanning_ft, 1)
            self.update_counts(reference_base, junction_wt1, spanning_wt1, 7)
            self.update_counts(reference_base, junction_wt2, spanning_wt2, 13)
        # set max anchor for this reference base
        self.fusion_seq_dict[reference_base][5] = anchor_ft
        self.fusion_seq_dict[reference_base][11] = anchor_wt1
        self.fusion_seq_dict[reference_base][17] = anchor_wt2
            #print("stored after: {}".format(self.fusion_seq_dict[reference_base]))

    def count_junc_span(self, read, breakpoint_pos, junction_count, spanning_counts, anchor):
        """Identify whether a read belongs to a junction or spanning pair"""
        # read overlaps the fusion breakpoint
        if breakpoint_pos >= self.bp_distance_threshold and read.get_overlap(breakpoint_pos - self.bp_distance_threshold, breakpoint_pos + self.bp_distance_threshold) == 2 * self.bp_distance_threshold:
            junction_count += 1
            anchor = max(anchor, min(read.reference_end - breakpoint_pos, breakpoint_pos - read.reference_start))
        # read maps left of the fusion breakpoint
        if read.reference_end < (breakpoint_pos + self.bp_distance_threshold):
            spanning_counts[0] += 1
        # read maps right of the fusion breakpoint
        if read.reference_start > (breakpoint_pos - self.bp_distance_threshold):
            spanning_counts[1] += 1
        return junction_count, spanning_counts, anchor

    # counter_start should be 1 for ft, 6 for wt1 and 11 for wt2 (see comment in header_to_dict for details)
    def update_counts(self, reference_base, junctions, spannings, counter_start):
        """Updates the read counter in the fusion_seq_dict"""
        if spannings[0] > 0:
            self.fusion_seq_dict[reference_base][counter_start] += 1
        if spannings[1] > 0:
            self.fusion_seq_dict[reference_base][counter_start + 1] += 1
        if junctions > 0:
            self.fusion_seq_dict[reference_base][counter_start + 2] += 1
        elif spannings[0] > 0 and spannings[0] == spannings[1]:
            self.fusion_seq_dict[reference_base][counter_start + 3] += 1

    def normalize_counts_cpm(self, count):
        """Returns the CPM of a read count based on the number of original input reads"""
        if self.input_read_count == 0:
            return count
        return count / float(self.input_read_count) * 1000000.0

    def header_to_dict(self):
        """Convert bam header into a dict with keys = SN fields and values = list of breakpoint + junction/spanning counter"""
        for header_seq_id in self.in_bam.header.to_dict()["SQ"]:
            # fusion dict counter
            # [0]: breakpoint position on ft
            # [1]: gene a counts on ft
            # [2]: gene b counts on ft
            # [3]: junction read count on ft
            # [4]: spanning pairs count on ft
            # [5]: longest anchor on ft
            
            # [6]: breakpoint position on wt1
            # [7]: gene a counts on wt1
            # [8]: gene b counts on wt2
            # [9]: junction read count on wt1
            # [10]: spanning pairs count on wt1
            # [11]: longest anchor on wt1
            
            # [12]: breakpoint position on wt2
            # [13]: gene a counts on wt2
            # [14]: gene b counts on wt2
            # [15]: junction read count on wt2
            # [16]: spanning pairs count on wt2
            # [17]: longest anchor on wt2
            
            # currently, the header is in the format: "FTID_contextSeqHash_breakpoint_{ft,wt1,wt2}"
            # As breakpoints can be different between ft and wt, basename is everything before
            id_split = header_seq_id["SN"].split("_")
            base_name = "_".join(id_split[:-2])
            # initialize the counter list for the respective fusion context
            if base_name not in self.fusion_seq_dict:
                self.fusion_seq_dict[base_name] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            if id_split[-1] == "ft":
                self.fusion_seq_dict[base_name][0] = int(id_split[-2])
            elif id_split[-1] == "wt1":
                self.fusion_seq_dict[base_name][6] = int(id_split[-2])
            elif id_split[-1] == "wt2":
                self.fusion_seq_dict[base_name][12] = int(id_split[-2])

        print("Found {} potential fusion sequences".format(len(self.fusion_seq_dict)))

    def run(self):
        """Walk linewise through a s/bam file and send reads mapping to the same fusion/wt context to counting"""
        #self.logger.info("Starting fusion read filtering")
        count_lines = 0
        count_processed_refs = 0
        # the read buffer is a dict of lists. The dict has query names as keys and a list of corresponding reads as values.
        read_buffer = {}
        last_reference = ""
        self.header_to_dict()

        for read in self.in_bam.fetch():
            count_lines += 1
            if not last_reference in read.reference_name and count_lines > 1:
                self.quantify_read_groups(read_buffer, last_reference)
                count_processed_refs += 1
                if count_processed_refs % (len(self.fusion_seq_dict) / 10) == 0:
                    print("At least {}% completed...".format(count_processed_refs / (len(self.fusion_seq_dict) / 100)))
                read_buffer.clear()
            if not read.query_name in read_buffer:
                read_buffer[read.query_name] = []
            read_buffer[read.query_name].append(read)
            # reference name conventions may change in the near future, but the tail must be
            # "*_breakpointPosition_ft" for the fusion transcript and "*_wt{1,2}" for the wt background
#            if read.reference_name.endswith("ft"):
            last_reference = "_".join(read.reference_name.split("_")[:-2])
#            else:
#                last_reference = "_".join(read.reference_name.split("_")[:-1])

        # process the last read_buffer
        self.quantify_read_groups(read_buffer, last_reference)
        print("All done! Reads did not map to ~{}% of the reference input".format(100 - (count_processed_refs / (len(self.fusion_seq_dict) / 100))))
        # close reading/writing stream
        self.in_bam.close()
        # write data to file. Change out_file_sep to ";" to create a csv file

        # get input read count
        with open(os.path.join(os.path.dirname(self.output), "Star_org_input_reads.txt"), "r") as rfile:
            self.input_read_count = int(rfile.next())

        out_file_sep = ";"
        header_string = out_file_sep.join(
            ["ftid_plus",
             "ft_bp", "ft_a", "ft_b", "ft_junc", "ft_span", "ft_anch",
             "wt1_bp", "wt1_a", "wt1_b", "wt1_junc", "wt1_span", "wt1_anch",
             "wt2_bp", "wt2_a", "wt2_b", "wt2_junc", "wt2_span", "wt2_anch"])
        # write counts and normalized counts
        # for normalization, breakpoint and anchor must not be converted: list positions 0, 5, 10, 15
        no_norm = {0, 5, 6, 11, 12, 17}
        with open("{}.counts".format(self.output), "w") as out_counts, open(self.output, "w") as out_norm:
            out_counts.write("{}\n".format(header_string))
            out_norm.write("{}\n".format(header_string))
            for key in self.fusion_seq_dict:
                out_counts.write("{}{}".format(key, out_file_sep))
                out_norm.write("{}{}".format(key, out_file_sep))
                # write counts directly
                out_counts.write("{}\n".format(out_file_sep.join(map(str, self.fusion_seq_dict[key]))))
                # normalize and write normalized counts
                self.fusion_seq_dict[key] = [read_count if i in no_norm else self.normalize_counts_cpm(read_count) for i, read_count in enumerate(self.fusion_seq_dict[key])]
                #self.fusion_seq_dict[key] = map(self.normalize_counts_cpm, self.fusion_seq_dict[key])
                out_norm.write("{}\n".format(out_file_sep.join(map(str, self.fusion_seq_dict[key]))))
        # write normalized counts
#        with
#            out_file2.write("{}\n".format(header_string))
#            for key in self.fusion_seq_dict:
                # write only those putative fusions, where not all counts are 0
                #if sum(self.fusion_seq_dict[key][1:]) > 0:
                #self.fusion_seq_dict[key].insert(0, key.rsplit("_", 2)[0])
#                self.fusion_seq_dict[key] = map(self.normalize_counts_cpm, self.fusion_seq_dict[key])
#                self.fusion_seq_dict[key].insert(0, key)
#                out_file.write("{}\n".format(out_file_sep.join(map(str, self.fusion_seq_dict[key]))))

#        out_file_sep = ";"
#        with open(self.output, "w") as out_file:
#            out_file.write(out_file_sep.join(["ftid_plus", "Type", "Reads_Gene_A", "Reads_Gene_B", "Junction_Reads", "Spanning_Pairs", "Longest_Anchor"]) + "\n")
#            for key in self.fusion_seq_dict:
#                # write only those putative fusions, where not all counts are 0
#                #if sum(self.fusion_seq_dict[key][1:]) > 0:
#                #ABCA7_19:1059141:+_ENST00000532194_DNHD1_11:6546666:+_ENST00000533649_9dd39c73a904776e_100_ft
#                ftid_p = key.rsplit("_", 2)[0]
#                (_, junc_ft, span_ft, anch_ft, junc_wt1, span_wt1, anch_wt1, junc_wt2, span_wt2, anch_wt2) = self.fusion_seq_dict[key]
#                out_file.write(out_file_sep.join(map(str, [ftid_p, "ft", "NA", "NA", junc_ft, span_ft, anch_ft])) + "\n")
#                out_file.write(out_file_sep.join(map(str, [ftid_p, "wt1", "NA", "NA", junc_wt1, span_wt1, anch_wt1])) + "\n")
#                out_file.write(out_file_sep.join(map(str, [ftid_p, "wt2", "NA", "NA", junc_wt2, span_wt2, anch_wt2])) + "\n")

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input', dest='input', help='Specify input BAM file', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Specify output file', required=True)
    parser.add_argument('-d', '--bp_distance', dest='bp_distance', help='Threshold of bases around the breakpoint for junction/spanning counting', default="3")
    args = parser.parse_args()

    requant = Requantification(args.input, args.output, args.bp_distance)
    requant.run()

if __name__ == '__main__':
    main()
