"""
Combining of all relevant infos into a final table
Calculations of additional metrics
Filtering of the final data table

@author: BNT (URLA)
@version: 20181121
"""

from __future__ import print_function
from argparse import ArgumentParser
import pandas as pd
import math
import glob
import os
import sys
from Bio import pairwise2
from easyfuse_iomethods import IOmethods as urlaio

# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class Fusionranking(object):
    """Select alignments belonging to putative fusions from an s/bam file"""
    def __init__(self, input_dir, id1, id2, output):
        """Parameter initialization"""        
        self.input_dir = input_dir
        self.id1 = id1
        self.id2 = id2
        self.output = output
        self.chars = ["A", "C", "G", "T"]
        pd.options.mode.chained_assignment = None

    @staticmethod
    def read_assembly_files(blast_assembly_dir):
        """Read de novo assembled fusion transcripts per gene fusion pair"""
        assembly_data = {}
        assembly_file_list = glob.glob(os.path.join(blast_assembly_dir, "*.blast.out.fusionTranscript"))
        for assembly_file in assembly_file_list:
            fusid = urlaio.path_leaf(assembly_file).split(".")[0]
            assembly_data[fusid] = {}
            with open(assembly_file, "r") as afile:
                tmp_header = ""
                for i, line in enumerate(afile, 1):
                    if i % 4 == 0:
                        tmp_header = line.split(" ")[0]
                        assembly_data[fusid][tmp_header] = ["NA", "NA", "NA"]
                    elif (i - 1) % 4 == 0:
                        assembly_data[fusid][tmp_header][0] = line.rstrip("\n")
                    elif (i - 2) % 4 == 0:
                        assembly_data[fusid][tmp_header][1] = line.rstrip("\n")
                    else:
                        assembly_data[fusid][tmp_header][2] = line.rstrip("\n")
        return assembly_data

    def entropy_calculation(self, sequence):
        """Return the shannon entropy of a sequence"""
        char_list = list(sequence)
        adj_seq_len = len(sequence)#-char_list.count("N")
        shannon_entropy = 0
        for c in self.chars:
            char_pc = float(char_list.count(c)) / float(adj_seq_len)
            if char_pc != 0:
                shannon_entropy += char_pc * math.log(char_pc, 2)
        return (-1) * shannon_entropy

    @staticmethod
    def similarity_calculation(sequence1, sequence2):
        """Return the score (relative to the max score) of a simple pairwise alignment with 1/-1/-2/-2 scores"""
        alignment = pairwise2.align.globalms(sequence1, sequence2, 1, -1, -2, -2)[0]
        return alignment[2] / float(max(len(sequence1), len(sequence2)))

    def filter_annotation_problems(self, pddf):
        """Identification of annotation problems probably as a result of different annotation versions"""
        count_put_miss_annot = 0
        for i in pddf.index:
            try:
                fusg1, fusg2 = pddf["Fusion_Gene"][i].split("_")
                if not fusg1 in pddf["FTID"][i] or not fusg2 in pddf["FTID"][i]:
                    count_put_miss_annot += 1
                    pddf = pddf.drop(i)
            except ValueError:
                count_put_miss_annot += 1
                pddf = pddf.drop(i)
        print("{} putative mis-annotations identified and dropped.".format(count_put_miss_annot))
        return pddf

    def append_to_context_seqs(self, pddf):
        """Append some metrics to the context seq data frame"""
        # urla: the dataframe will require further filtering to be implemented which will be easy starting with a pd df.
        #       furthermore, it garantues to generate a valid csv file.
        #win test:C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\innoData\\fusionProcessing_optV2fltr\\Context_Seqs.csv
        # C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\PRJNA252360_data\\100bp_filtered\\SRR1659951\\classification.tdt
        # C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\PRJNA252360_data\\100bp_filtered\\SRR1659951\\Context_Seqs.csv
        pddf["left_bp_entropy"] = ""
        pddf["right_bp_entropy"] = ""
        pddf["bp_homology_score"] = ""
        for i in pddf.index:
            bp1_sequence = pddf["context_sequence"][i][pddf["context_sequence_bp"][i]:]
            bp2_sequence = pddf["context_sequence"][i][:pddf["context_sequence_bp"][i]]
            pddf["left_bp_entropy"][i] = self.entropy_calculation(bp1_sequence)
            pddf["right_bp_entropy"][i] = self.entropy_calculation(bp2_sequence)
            pddf["bp_homology_score"] = self.similarity_calculation(bp1_sequence, bp2_sequence)
        return pddf
        #pdcsv.to_csv(self.output, sep = "\t", index=False)

    def read_cnt_filtering(self, pddf):
        """Keep rows fulfilling the defined criteria and drop the rest"""
        for i in pddf.index:
            count_checks = [0, 0, 0, 0, 0, 0, 0, 0]
            # if fusion is supported in sample 1 or 2
            if pddf["junction_ft_1"][i] > 0 and pddf["spanning_ft_1"][i] > 0:
                count_checks[0] = 1
            if pddf["junction_ft_2"][i] > 0 and pddf["spanning_ft_2"][i] > 0:
                count_checks[0] *= 3
                count_checks[1] = 1
            
            # if wt background is not supported in sample 1 or 2
            if pddf["junction_wt_1"][i] == 0 and pddf["spanning_wt_1"][i] == 0:
                count_checks[2] = 1
            if pddf["junction_wt_2"][i] == 0 and pddf["spanning_wt_2"][i] == 0:
                count_checks[2] *= 3
                count_checks[3] = 1
                
            # if fusion is support by a long (>= 20) anchor read while wt is not
            if pddf["longest_anchor_ft_1"][i] >= 20 and pddf["longest_anchor_wt_1"][i] < 20:
                count_checks[4] = 1
            if pddf["longest_anchor_ft_2"][i] >= 20 and pddf["longest_anchor_wt_2"][i] < 20:
                count_checks[5] = 1
            
            # if fusion is on exon boundary
            if pddf["exon_boundary"][i] == "both":
                count_checks[6] = 2
            elif pddf["exon_boundary"][i] == "none":
                count_checks[6] = 0
            else:
                count_checks[6] = -2

            # if fusion transcript is out of frame
            if pddf["frame"][i] == "no_frame":
                count_checks[7] = -99
            else:
                count_checks[7] = 2

            if sum(count_checks) < 10:
                pddf = pddf.drop(i)

        return pddf

    def stability_calculation(self, sequence, dipeptide_score_matrix):
        """Retrun the stability index as proposed by Guruprasad et al., 1990"""
        dsm = pd.read_table(dipeptide_score_matrix, header=0, index_col=0)
        stab_idx = 0
        for i in range(len(sequence) - 1):
            stab_idx += dsm[sequence[i+1]][sequence[i]]
        return stab_idx * 10 / float(len(sequence))

    def run_one(self):
        print("Looking for required files...")
        #win test
        context_seq_file1 = "C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\icamData_immunTest\\A12656_R1_Context_Seqs.csv"
        context_seq_file2 = "C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\icamData_immunTest\\A12656_R2_Context_Seqs.csv"
        requant_file1 = "C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\icamData_immunTest\\A12656_R1_classification.tdt"
        requant_file2 = "C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\icamData_immunTest\\A12656_R2_classification.tdt"
        # context seq files
        context_seq_file1 = os.path.join(self.input_dir, "Sample_{}".format(self.id1), "fetchdata", "fd_1_tool", "fetched_contextseqs", "Context_Seqs.csv")
        self.check_files(context_seq_file1)
        context_seq_file2 = os.path.join(self.input_dir, "Sample_{}".format(self.id2), "fetchdata", "fd_1_tool", "fetched_contextseqs", "Context_Seqs.csv")
        self.check_files(context_seq_file2)
        # re-quantification files
        requant_file1 = os.path.join(self.input_dir, "Sample_{}".format(self.id1), "fetchdata", "fd_1_tool", "classification", "classification.tdt")
        self.check_files(requant_file1)
        requant_file2 = os.path.join(self.input_dir, "Sample_{}".format(self.id2), "fetchdata", "fd_1_tool", "classification", "classification.tdt")
        self.check_files(requant_file2)
        # assembly data
        assembly_dir1 = os.path.join(self.input_dir, "Sample_{}".format(self.id1), "fetchdata", "fd_1_tool", "assembly")
        #self.check_files(assembly_dir1)
        assembly_dir2 = os.path.join(self.input_dir, "Sample_{}".format(self.id2), "fetchdata", "fd_1_tool", "assembly")
        #self.check_files(assembly_dir2)

        print("Loading data into memory...")
        context_data1 = pd.read_csv(context_seq_file1, sep = ";")
        context_data1 = self.filter_annotation_problems(context_data1)
        context_data1["ftid_plus"] = context_data1["FTID"] + "_" +  context_data1["context_sequence_id"] # appending key
        requant_data1 = pd.read_table(requant_file1)
        requant_data2 = pd.read_table(requant_file2)
        assembl_data1 = self.read_assembly_files(assembly_dir1)
        assembl_data2 = self.read_assembly_files(assembly_dir2)

        print("Merging re-quantification data (inner join)")
        # urla: I start from the re-quantification because it provides the most meaningful filter and is the easiest data to process
        requant_data12 = requant_data1.set_index('fusion_id').join(requant_data2.set_index('fusion_id'), lsuffix = "_1", rsuffix = "_2", how = 'inner')
#        requant_data12.to_csv(self.output, sep = "\t", index=True)

        print("Merging re-quantification data into context seq (inner join)")
        context_data12 = context_data1.set_index('ftid_plus').join(requant_data12, how = 'inner')
        context_data12.to_csv(self.output, sep = "\t", index=True)
        
        print("Appending entropy and homology calculations")
        context_data12 = self.append_to_context_seqs(context_data12)
        context_data12.to_csv(self.output, sep = "\t", index=True)

        print("Performing read support filtering and preparing de novo assembly targets")
        context_data12_filtered = self.read_cnt_filtering(context_data12)
        print("Cleaning final data")
        context_data12_filtered = context_data12_filtered.drop(
                ['FGID',
                 'Fusion_Gene',
                 'Breakpoint1',
                 'Breakpoint2',
                 'FTID',
                 'context_sequence_id',
                 'exon_boundary1',
                 'exon_boundary2',
                 'wt1_context_sequence',
                 'wt1_context_sequence_bp',
                 'wt2_context_sequence',
                 'wt2_context_sequence_bp',
                 'breakpoint_pos_1',
                 'longest_anchor_ft_1',
                 'longest_anchor_wt_1',
                 'breakpoint_pos_2',
                 'longest_anchor_ft_2',
                 'longest_anchor_wt_2'
                 ], axis = 1
                 )
        context_data12_filtered.to_csv(self.output + ".filtered", sep = "\t", index=True)
        
#        with open(os.path.join(assembly_dir1, "fusion_gene_list.txt"), "w") as fgl1, open(os.path.join(assembly_dir2, "fusion_gene_list.txt"), "w") as fgl2:
#            for fusion_pair in set(context_data12["Fusion_Gene"]):
#                fgl1.write("{}\n".format(fusion_pair.replace("_", "--")))
#                fgl2.write("{}\n".format(fusion_pair.replace("_", "--")))
        

    def run_two(self):
        # assembly data
        assembly_dir1 = os.path.join(self.input_dir, "Sample_{}".format(self.id1), "fetchdata", "fd_1_tool", "assembly", "blast")
        #self.check_files(assembly_dir1)
        assembly_dir2 = os.path.join(self.input_dir, "Sample_{}".format(self.id2), "fetchdata", "fd_1_tool", "assembly", "blast")
        #self.check_files(assembly_dir2)
        #context_data12 = pd.read_table(self.output)
        print("Reading de novo assembly data")

    @staticmethod
    def check_files(file_path):
        if not os.path.isfile(file_path) and not os.path.isdir(file_path):
            print("File or directory not found or just not at the expected location. Looked at {}".format(file_path))
            sys.exit(1)
        else:
            print("Found {}".format(file_path))

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('-i', '--input_dir', dest='input_dir', help='Input root directory containing the sample file directories', required=True)
    parser.add_argument('--id1', dest='id1', help='Identifier of the first sample replicate', required=True)
    parser.add_argument('--id2', dest='id2', help='Identifier of the second sample replicate', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Specify output file', required=True)
    parser.add_argument('--do2not1', dest='do2not1', help='Run second filtering (assembly)', default=False, action='store_true')
    args = parser.parse_args()

    fusrank = Fusionranking(args.input_dir, args.id1, args.id2, args.output)
    if args.do2not1:
        fusrank.run_two()
    else:
        fusrank.run_one()

if __name__ == '__main__':
    main()
