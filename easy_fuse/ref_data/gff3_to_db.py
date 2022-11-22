#!/usr/bin/env python

"""
Create (if necessary) a sqlite3 database from an ensembl gff3 file.
Afterwards, test whether the db can be accessed and print the counts
of all features in it

@author: BNT (URLA)
@version: 20190704
"""

from __future__ import print_function, division
import os
import argparse
import gffutils


# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class GTF2DB(object):
    """stub"""
    def __init__(self, gff_in_file, db_out_file):
        """Parameter initialization"""
        self.gff_in_file = gff_in_file
        self.db_out_file = db_out_file

    def create_gffdb(self, overwrite):
        """ Create a gffutils database from an ensembl gff3 file """
        if os.path.exists(self.db_out_file):
            if overwrite:
                print("File {} already exists, but overwrite is specified. Trying to delete the existing file".format(self.db_out_file))
                os.remove(self.db_out_file)
                self.create_gffdb(overwrite)
            else:
                print("File {} already exists and overwrite is NOT specified.".format(self.db_out_file))
        else:
            print("File {} does not exist. Trying to create it from the provided gff file (This may take 10-20+ minutes)...".format(self.db_out_file))
            # this is a one-time operation and tailored to human ensembl gff3
            gffutils.create_db(self.gff_in_file, self.db_out_file, force_gff=True, force=True, merge_strategy="create_unique")

        print("Testing whether db can be accessed...")
        self.test_db()

    def test_db(self):
        """ Test whether an existing database can be loaded and accessed """
        gff3_db = gffutils.FeatureDB(self.db_out_file)
        for feature in gff3_db.featuretypes():
            print("Found {} features of type {} in the db.".format(gff3_db.count_features_of_type(feature), feature))

def main():
    """Parse command line arguments and start script"""
    parser = argparse.ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument('--gff_input', dest='gff_in', help='ensembl gff3 annotation file', required=True)
    parser.add_argument('--db_output', dest='db_out', help='gff3 derived annotation database from gffutils', required=True)
    parser.add_argument('--overwrite', dest='overwrite', help='Overwrite existing database if present', default=False, action='store_true')
    args = parser.parse_args()
    GTF2DB(args.gff_in, args.db_out).create_gffdb(args.overwrite)

if __name__ == '__main__':
    main()
