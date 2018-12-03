"""
Junction/Spanning read pair counting

@author: BNT (URLA)
@version: 20181126
"""

from __future__ import print_function
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

    def quantify_read_groups(self, read_buffer, reference_base):
        """Count junction/spanning pairs on ft and wt1/2"""
        breakpoint_pos = self.fusion_seq_dict[reference_base][0]
        anchor_ft = 0
        anchor_wt = 0
        for read_group in read_buffer:
            junction_ft = 0 # 0 = no junction overlapping read, 1+ = at least one of the paired reads is a junction read
            spanning_ft = [0, 0] # left: paired read cnt fully left of the fusion bp, right: paired read cnt fulls right of the fusion bp
            junction_wt = 0 # like is_junction_ft but for wt
            spanning_wt = [0, 0] # like is_spanning_ft but for wt

            for read in read_buffer[read_group]:
                # exclude non-primary alignments
                #if read.flag > 255:
                #    continue
                # exclude low mapq alignments
                #if read.mapping_quality < 1:
                #    continue
                if read.reference_name.endswith("ft"):
                    junction_ft, spanning_ft, anchor_ft = self.count_junc_span(read, breakpoint_pos, junction_ft, spanning_ft, anchor_ft)
                else:
                    junction_wt, spanning_wt, anchor_wt = self.count_junc_span(read, breakpoint_pos, junction_wt, spanning_wt, anchor_wt)

            #print("stored before: {}".format(self.fusion_seq_dict[reference_base]))
            #print("jft: {}, sft: {}, jwt: {}, swt: {}".format(is_junction_ft, is_spanning_ft, is_junction_wt, is_spanning_wt))
            if junction_ft > 0:
                self.fusion_seq_dict[reference_base][1] += 1
            elif spanning_ft[0] > 0 and spanning_ft[0] == spanning_ft[1]:
                self.fusion_seq_dict[reference_base][2] += 1
            if junction_wt > 0:
                self.fusion_seq_dict[reference_base][3] += 1
            elif spanning_wt[0] > 0 and spanning_wt[0] == spanning_wt[1]:
                self.fusion_seq_dict[reference_base][4] += 1
        self.fusion_seq_dict[reference_base][5] = anchor_ft
        self.fusion_seq_dict[reference_base][6] = anchor_wt
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

    def header_to_dict(self):
        """Convert bam header into a dict with keys = SN fields and values = list of breakpoint + junction/spanning counter"""
        for header_seq_id in self.in_bam.header.to_dict()["SQ"]:
            if header_seq_id["SN"].endswith("_ft"):
                id_split = header_seq_id["SN"].split("_")
                self.fusion_seq_dict["_".join(id_split[:-2])] = [int(id_split[-2]), 0, 0, 0, 0, 0, 0] # initialize with bp loc on ft as well as counter for junction ft, spanning ft, junction wt, spanning wt mappings and longest anchor
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
            if read.reference_name.endswith("ft"):
                last_reference = "_".join(read.reference_name.split("_")[:-2])
            else:
                last_reference = "_".join(read.reference_name.split("_")[:-1])

        # process the last read_buffer
        self.quantify_read_groups(read_buffer, last_reference)
        print("All done! Reads did not map to ~{}% of the reference input".format(100 - (count_processed_refs / (len(self.fusion_seq_dict) / 100))))
        # close reading/writing stream
        self.in_bam.close()
        # write data to file. Change out_file_sep to ";" to create a csv file
        out_file_sep = "\t"
        with open(self.output, "w") as out_file:
            out_file.write(out_file_sep.join(["fusion_id", "breakpoint_pos", "junction_ft", "spanning_ft", "junction_wt", "spanning_wt", "longest_anchor_ft", "longest_anchor_wt"]) + "\n")
            for key in self.fusion_seq_dict:
                # write only those putative fusions, where not all counts are 0
                #if sum(self.fusion_seq_dict[key][1:]) > 0:
                self.fusion_seq_dict[key].insert(0, key)
                out_file.write(out_file_sep.join(map(str, self.fusion_seq_dict[key])) + "\n")

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
