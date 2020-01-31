#!/usr/bin/env python

"""
Sample file handling
Sample processing is logged in a sqlite3 database. This class provides methods to
create and work on this DB. Concerning data structures: All data is stored
in a samples table including sample_id, list of performed analyses as well as the corresponding FASTQ file location.
Data itself can be modified using wrapped SQL queries.

<<<<<<< HEAD
@author: Tron (PASO), BNT (URLA)
@version: 20190429
"""


from argparse import ArgumentParser

import sqlite3

class SamplesDB(object):
    """Create and interact with the sample DB"""
    def __init__(self, db_path):
        """Parameter initialization and connection to database"""
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.cursor = self.conn.cursor()
        self.create_table()

    def create_table(self):
        """Creation of samples table if non-existent"""
        self.cursor.execute("""CREATE TABLE IF NOT EXISTS samples (sample_id text, analysis text, fq1 text, fq2 text)""")
        self.commit_changes()

    def add_sample(self, sample_id, analysis, fq1, fq2):
        """Add sample to DB if not exists"""
        if self.get_sample(sample_id):
            return
        else:
            self.cursor.execute("INSERT INTO samples VALUES ('{}','{}','{}','{}')".format(sample_id, analysis, fq1, fq2))
            self.commit_changes()

    def get_sample(self, sample_id):
        """Get sample entry from DB"""
        self.cursor.execute("SELECT * FROM samples WHERE sample_id = '{}'".format(sample_id))
        sample_id = self.cursor.fetchone()
        if sample_id:
            return sample_id[0]
        else:
            return None
=======
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
>>>>>>> 55fa52a168eebc0ae4f190e4019c00a654ba9fa0

    def get_tool_list_from_state(self, sample_id):
        """Get list of performed analyses from sample"""
        self.cursor.execute("SELECT analysis FROM samples WHERE sample_id = '{}'".format(sample_id))
        tool_list = self.cursor.fetchone()[0].split(",")
        return tool_list

    def get_fastq_files(self, sample_id):
<<<<<<< HEAD
        """Get fastq files from sample"""
        self.cursor.execute("SELECT fq1, fq2 FROM samples WHERE sample_id = '{}'".format(sample_id))
        (fq1, fq2) = self.cursor.fetchone()
        return (fq1, fq2)

    def get_sample_ids(self):
        """Get all sample_ids from DB"""
        sample_ids = []
        for row in self.cursor.execute("SELECT sample_id FROM samples"):
            sample_ids.append(row[0])
        return sample_ids

    def get_samples(self):
        """Get all existing samples in DB"""
        self.cursor.execute("SELECT * FROM samples")
        return self.cursor.fetchall()

    def print_db(self):
        """Print sample progress"""
        completed = 0
        total = 0
        for row in self.cursor.execute("SELECT * FROM samples"):
            print("Sample ID:", str(row[0]))
            print("Performed Analyses:", str(row[1]))
            print("FQ1:", str(row[2]))
            print("FQ2:", str(row[3]))
            if "Fetchdata" in row[1]:
                print("Status: Complete")
                completed += 1
            else:
                print("Status: Incomplete")
            total += 1
            print
        print("Overall: {:.2f}% complete.".format(completed*100.0/total))

    def append_state(self, sample_id, analysis):
        """Add analysis to list of performed analyses within sample"""
        curr_analysis = self.get_tool_list_from_state(sample_id)
        if analysis in curr_analysis:
            return
        if curr_analysis[0] == "NA":
            analyses = analysis
        else:
            analyses = "{},{}".format(",".join(curr_analysis), analysis)
        self.cursor.execute("UPDATE samples SET analysis = '{}' WHERE sample_id = '{}'".format(analyses, sample_id))
        self.commit_changes()

    def commit_changes(self):
        """Commit changes to DB connection"""
        self.conn.commit()

    def close_connection(self):
        """Close DB connection"""
        self.conn.close()
=======
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
>>>>>>> 55fa52a168eebc0ae4f190e4019c00a654ba9fa0

def main():
    """Main method, parameter handling"""
    parser = ArgumentParser(description='Handle sample db operations')

    parser.add_argument('-i', '--sample_id', dest='sample_id', help='Specify the sample id to process.', required=True)
    parser.add_argument('-d', '--db_path', dest='db_path', help='Specify the file to save the changes into.', required=True)
    parser.add_argument('-t', '--tool', dest='tool', help='The tools name to add to your sample id.', required=True)
    parser.add_argument('-l', '--action', dest='action', choices=['append_state'], help='Select the action to do with the sample id.', required=True)
    args = parser.parse_args()
    
    db = SamplesDB(args.db_path)

    if args.action == "append_state":
<<<<<<< HEAD
        db.append_state(args.sample_id, args.tool)
    db.close_connection()
=======
        samples.append_state(args.sample_id, args.state)
    elif args.action == "print_db":
        samples.write_file()
    elif args.action == "convert_db":
        samples.csv2db()
>>>>>>> 55fa52a168eebc0ae4f190e4019c00a654ba9fa0

if __name__ == "__main__":
    main()
