#!/usr/bin/env python

"""
Sample file handling
Sample processing is logged in a specific file. This class provides methods to
create and work on this file. Concerning data structures: All data is stored
in a dict with sample identifiers as keys and data as value. Data itself is
are 2-tuples of lists, currently including an identifier of the tools that have been
run on the sample and the fastq files corresponding to the sample.

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""


from __future__ import print_function
import os
from argparse import ArgumentParser

class Samples(object):
    """Create and interact with the sample monitoring file"""
    def __init__(self, infile):
        """Parameter initialization"""
        self.infile = infile
        if not os.path.exists(self.infile):
            open(self.infile, "a").close()
        self.sample_map = self.read_file()

    # getter and setter methods
    def get_state(self, sample_id):
        """Get the state of a sample based on its id"""
        state = ""
        if sample_id in self.sample_map:
            state = self.sample_map[sample_id][0]
        return state

    def set_state(self, sample_id, state):
        """Update, that is override, a sample plus corresponding state to a sample file"""
        if sample_id in self.sample_map:
            fastq = self.sample_map[sample_id][1]
            self.sample_map[sample_id] = (state, fastq)
        self.write_file()

    def get_fastq_files(self, sample_id):
        """Get a tuple of fastq file strings"""
        fastq = ""
        if sample_id in self.sample_map:
            fastq = self.sample_map[sample_id][1]
            if "R2" in fastq:
                return (fastq.split(";")[0][3:], fastq.split(";")[1][3:])
            return fastq[3:]
        return fastq

    def set_fastq_files(self, sample_id, fastq):
        """Set a tuple of fastq file strings"""
        if sample_id in self.sample_map:
            state = self.sample_map[sample_id][0]
            self.sample_map[sample_id] = (state, fastq)
        self.write_file()


    def add_sample(self, sample_id, state, fastq):
        """Add a sample plus corresponding state and fastq file(s) to a sample file"""
        if not sample_id in self.sample_map:
            self.sample_map[sample_id] = (state, fastq)
        self.write_file()

    def delete_sample(self, sample_id):
        """Delete a sample"""
        if sample_id in self.sample_map:
            del self.sample_map[sample_id]
        self.write_file()

    def append_state(self, sample_id, state):
        """Append something to an existing state or use update"""
        if sample_id in self.sample_map:
            if self.get_state_string(sample_id) == "NA":
                self.set_state(sample_id, state)
            elif state not in self.get_state_string(sample_id):
                self.sample_map[sample_id] = ("{0}, {1}".format(self.get_state(sample_id), state), self.sample_map[sample_id][1])
        self.write_file()


    # read and write sample file
    def read_file(self):
        """Read a sample file and return a dictionary (sample_id: (state, fastq files))"""
        sample_map = {}
        with open(self.infile) as sample_file:
            next(sample_file)
            for line in sample_file:
                elements = line.rstrip().split("\t")
                try:
                    sample_id = elements[0]
                    state = elements[1]
                    fastq = elements[2]
                    sample_map[sample_id] = (state, fastq)
                except IndexError as idx_err:
                    print(idx_err)
                    print(line)
        return sample_map

    def write_file(self):
        """Write a sample file to disk"""
        with open(self.infile, 'w') as outf:
            outf.write("{}\t{}\n".format("#sample_id", "state", "fastq_files"))
            for key in self.sample_map:
                (state, fastq) = self.sample_map[key]
                outf.write("{}\t{}\t{}\n".format(key, state, fastq))


def main():
    """Main method, parameter handling"""
    parser = ArgumentParser(description='Handle sample db operations')
                       
    parser.add_argument('-i', '--sample-id', dest='sample_id', help='Specify the sample id to process.', required=True)
    parser.add_argument('-o', '--output-file', dest='output_file', help='Specify the file to save the changes into.', required=True)
    parser.add_argument('-g', '--state', dest='state', help='Specify the new state of your sample id.', required=True)
    parser.add_argument('-l', '--action', dest='action', choices=['append_state'], help='Select the action to do with the sample id.', required=True)
    args = parser.parse_args()

    samples = Samples(args.output_file)
    if args.action == "append_state":
        samples.append_state(args.sample_id, args.state)
                       


if __name__ == '__main__':
    main()
