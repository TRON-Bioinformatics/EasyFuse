#!/usr/bin/env python

"""
Sample file handling
Sample processing is logged in a specific file. This class provides methods to
create and work on this file. Concerning data structures: All data is stored
in a dict with sample identifiers as keys and data as value. Data itself is
are 2-tuples of lists, currently including an identifier of the tools that have been
run on the sample and the fastq files corresponding to the sample.

@author: BNT (URLA)
@version: 20190507
"""

from __future__ import print_function
import sqlite3
from argparse import ArgumentParser

# pylint: disable=line-too-long
#         urla: yes they are partially, but I do not consider this to be relevant here
class Samples(object):
    """Create and interact with the sample monitoring file"""
    def __init__(self, infile):
        """Parameter initialization"""
        self.infile = infile
        self.db_con = sqlite3.connect(self.infile)
        # check if the table "samples" exists and create it if not
        try:
            self.db_con.execute("create table samples (sample_id varchar primary key, state varchar, fastq_file varchar)")
        except sqlite3.OperationalError:
            # table already exists
            pass

    # getter and setter methods
    def get_sample_id_list(self):
        """ Return a list of sample_ids from the database """
        db_cursor = self.db_con.cursor()
        ids = db_cursor.execute("SELECT sample_id FROM samples").fetchall()
        return [sid for (sid,) in ids]

    def get_state_string(self, sample_id):
        """Get the state of a sample based on its id"""
        db_cursor = self.db_con.cursor()
        try:
            return db_cursor.execute("SELECT state FROM samples WHERE sample_id='{}';".format(sample_id)).fetchone()[0]
        except TypeError:
            print("Sample {} is not available in the database".format(sample_id))
            return ""

    def get_tool_list_from_state(self, sample_id):
        """Get a list of tools that have been successfully run on a sample based on its id"""
        state = self.get_state_string(sample_id)
        tools = []
        if ", " in state:
            for tool in state.split(", "):
                tools.append(tool.split("-")[0])
        else:
            tools = state.split("-")[0]
        return tools

    def get_fastq_files(self, sample_id):
        """Get a tuple of fastq file strings"""
        fastq = ""
        db_cursor = self.db_con.cursor()
        try:
            fastq = db_cursor.execute("SELECT fastq_file FROM samples WHERE sample_id='{}';".format(sample_id)).fetchone()[0]
            if "R2" in fastq:
                return (fastq.split(";")[0][3:], fastq.split(";")[1][3:])
            return fastq[3:]
        except TypeError:
            print("Sample {} is not available in the database".format(sample_id))
            return fastq

    def add_sample(self, sample_id, state, fastq):
        """Add a sample plus corresponding state and fastq file(s) to a sample file"""
        try:
            with self.db_con:
                self.db_con.execute("INSERT INTO samples(sample_id, state, fastq_file) values (?,?,?)", (sample_id, state, fastq))
        except sqlite3.IntegrityError:
            print("Sample {} is already in the table".format(sample_id))

    def set_state(self, sample_id, state):
        """Update, that is override, a sample plus corresponding state to a sample file"""
        with self.db_con:
            self.db_con.execute("UPDATE samples SET state = ? WHERE sample_id = ?", (state, sample_id))

    def append_state(self, sample_id, state):
        """Append something to an existing state or use update"""
        state_str = self.get_state_string(sample_id)
        if state_str == 'NA':
            self.set_state(sample_id, state)
        elif state not in state_str:
            self.set_state(sample_id, "{0}, {1}".format(state_str, state))

    def write_file(self):
        """ Write the whole db to file """
        db_cursor = self.db_con.cursor()
        db_out = db_cursor.execute("SELECT * FROM samples;").fetchall()

        out_file = "{}.csv".format(self.infile[:-3])
        with open(out_file, 'w') as out:
            out.write("#Sample ID\tPerformed analyses\tCorresponding fastq files\n")
            for line in ["\t".join([x, y, z]) for (x, y, z) in db_out]:
                out.write("{}\n".format(line))

    def csv2db(self):
        """ This method allows conversion from formerly used samples.csv to sqlite3 samples.sb files """
        in_csv = "{}.csv".format(self.infile[:-3])
        try:
            with self.db_con:
                with open(in_csv, "r") as in_file:
                    for line in in_file:
                        if not line or line.startswith("#"):
                            continue
                        line_splitter = line.rstrip().split("\t")
                        self.db_con.execute("INSERT INTO samples(sample_id, state, fastq_file) values (?,?,?)", (line_splitter[0], line_splitter[1], line_splitter[2]))
        except sqlite3.IntegrityError:
            print("Sample {} is already in the db".format(line_splitter[0]))

def main():
    """Main method, parameter handling"""
    parser = ArgumentParser(description='Handle sample db operations')

    parser.add_argument('-i', '--sample-id', dest='sample_id', help='Specify the sample id to process.', required=False)
    parser.add_argument('-o', '--output-file', dest='output_file', help='Specify the file to save the changes into.', required=True)
    parser.add_argument('-g', '--state', dest='state', help='Specify the new state of your sample id.', required=False)
    parser.add_argument('-l', '--action', dest='action', choices=['append_state', 'print_db', 'convert_db'], help='Select the action to do with the sample id.', required=True)
    args = parser.parse_args()

    samples = Samples(args.output_file)
    if args.action == "append_state":
        samples.append_state(args.sample_id, args.state)
    elif args.action == "print_db":
        samples.write_file()
    elif args.action == "convert_db":
        samples.csv2db()

if __name__ == '__main__':
    main()
