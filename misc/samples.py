#!/usr/bin/env python

"""
Sample file handling
Sample processing is logged in a sqlite3 database. This class provides methods to
create and work on this DB. Concerning data structures: All data is stored
in a samples table including sample_id, list of performed analyses as well as the corresponding FASTQ file location.
Data itself can be modified using wrapped SQL queries.

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

    def get_tool_list_from_state(self, sample_id):
        """Get list of performed analyses from sample"""
        self.cursor.execute("SELECT analysis FROM samples WHERE sample_id = '{}'".format(sample_id))
        tool_list = self.cursor.fetchone()[0].split(",")
        return tool_list

    def get_fastq_files(self, sample_id):
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
        db.append_state(args.sample_id, args.tool)
    db.close_connection()

if __name__ == "__main__":
    main()
