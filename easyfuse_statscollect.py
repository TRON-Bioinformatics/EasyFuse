"""
Sample stats handling
Stats are organised in a dict of dict of dicts.
outer dict: sample ids
middle dict: tool ids
inner dict: collected stats id
value: collected stats value
The aim of this method is to preserver stats infos which may not be saved elsewhere
(e.g. runtimes provided by sacct)

@author: BNT (URLA)
@version: 20181126
"""
import os
from argparse import ArgumentParser

class SCollector(object):
    """Create and interact with the sample monitoring file"""
    def __init__(self, infile):
        """Parameter initialization"""
        self.infile = infile
        if not os.path.exists(self.infile):
            open(self.infile, "a").close()
        self.data_map = self.read_file()

    # read and write sample file
    def read_file(self):
        """Read a sample file and return a dictionary (sample_id; (state, fastq files))"""
        sample_map = {}
        tools_map = {}
        stats_map = {}
        with open(self.infile) as sample_file:
            for i, line in enumerate(sample_file):
         #       print("l{}: {}".format(i, line))
                line_splitter = line.rstrip().split("\t")
          #      print("ls: {}".format(line_splitter))
                if len(line_splitter) > 0 and not line.startswith("#"):
                    if not line_splitter[0] == "":
                        sample_id = line_splitter[0]
                        tools_map = {}
                        stats_map = {}
                    if not line_splitter[1] == "":
                        tool_uid = line_splitter[1]
                        stats_map = {}
                    stats_id = line_splitter[2]
                    stats_val = line_splitter[3]
                    stats_map[stats_id] = stats_val
                    tools_map[tool_uid] = stats_map
                    sample_map[sample_id] = tools_map
        #print(sample_map)
        return sample_map

    def write_file(self):
        """Write a sample file to disk"""
        #print(self.data_map)
        with open(self.infile, "w") as outf:
            outf.write("{0}\t{1}\t{2}\t{3}\n".format("#Sample ID", "Tool UID", "Collected Stat", "Stats val"))
            for sample_id in self.data_map:
                for i, tool_uid in enumerate(self.data_map[sample_id]):
                    for j, stats_id in enumerate(self.data_map[sample_id][tool_uid]):
                        if i == 0 and j == 0:
                            outf.write("{0}\t{1}\t{2}\t{3}\n".format(sample_id, tool_uid, stats_id, self.data_map[sample_id][tool_uid][stats_id]))
                        elif j == 0:
                            outf.write("{0}\t{1}\t{2}\t{3}\n".format("", tool_uid, stats_id, self.data_map[sample_id][tool_uid][stats_id]))
                        else:
                            outf.write("{0}\t{1}\t{2}\t{3}\n".format("", "", stats_id, self.data_map[sample_id][tool_uid][stats_id]))

    # getter and setter methods
    def get_stat_from_tool(self, tool_id, stats_id):
        """Get a list of tools that have been successfully run on a sample based on its id"""
        stats = []
        for sample_id in self.data_map:
            for tool_uid in self.data_map[sample_id]:
                if tool_id in tool_uid:
                    if stats_id in self.data_map[sample_id][tool_uid]:
                        stats.append((sample_id, tool_uid, stats_id, self.data_map[sample_id][tool_uid][stats_id]))
        return stats

    def add_sample(self, sample_id):
        """Add a sample plus corresponding state and fastq file(s) to a sample file"""
        if not sample_id in self.data_map:
            self.data_map[sample_id] = {}

    def add_tool(self, sample_id, tool_uid):
        """Add a sample plus corresponding state and fastq file(s) to a sample file"""
        self.add_sample(sample_id)
        if not tool_uid in self.data_map[sample_id]:
            self.data_map[sample_id][tool_uid] = {}

    def add_stat(self, sample_id, tool_uid, stats_id, stats_val):
        """Append something to an existing state or use update"""
        self.add_tool(sample_id, tool_uid)
        self.data_map[sample_id][tool_uid][stats_id] = stats_val
        self.write_file()

def main():
    """Main method, parameter handling"""
    parser = ArgumentParser(description='Handle sample db operations')

    parser.add_argument('-i', '--sample-id', dest='sample_id', help='Specify the sample id to process.')
    parser.add_argument('-t', '--tool-id', dest='tool_id', help='Specify the sample id to process.')
    parser.add_argument('-o', '--output-file', dest='output_file', help='Specify the file to save the changes into.', required=True)
    parser.add_argument('-s', '--stat-id', dest='stat_id', help='Specify the new state of your sample id.')
    parser.add_argument('-a', '--add-val', dest='stat_val', help='Specify the new state of your sample id.')
    parser.add_argument('-g', '--get', dest='get_val', help='Specify the new state of your sample id.')
    args = parser.parse_args()

    stats = SCollector(args.output_file)
    #stats.add_stat(args.sample_id, args.tool_id, args.stat_id, args.stat_val)
    print(stats.get_stat_from_tool(args.tool_id, args.stat_id))

if __name__ == '__main__':
    main()
