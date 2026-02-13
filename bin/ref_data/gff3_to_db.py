#!/usr/bin/env python

"""
Create (if necessary) a sqlite3 database from an ensembl gff3 file.
Afterwards, test whether the db can be accessed and print the counts
of all features in it

@author: TRON (PASO)
@version: 20241119
"""

from argparse import ArgumentParser
import logging
import os

# pylint: disable=E0401
import gffutils # type: ignore

# FORMAT = '%(asctime)s - %(message)s'
# logging.basicConfig(format=FORMAT, level=logging.INFO)
# logger = logging.getLogger(__name__)


class Gff2Db:
    """Class for creating a gffutils database from a gff3 file"""
    def __init__(self, gff_in_file: str, db_out_file: str):
        """Parameter initialization"""
        self.gff_in_file = gff_in_file
        self.db_out_file = db_out_file


    def clean_overwrite(self, overwrite: bool):
        """Delete existing database if overwrite is specified"""
        if os.path.exists(self.db_out_file) and overwrite:
            # logger.info("File %s already exists, but overwrite is specified.", self.db_out_file)
            # logger.info("Deleting file %s", self.db_out_file)
            os.remove(self.db_out_file)
        else:
            pass
            # logger.info("File %s already exists or overwrite is NOT specified.", self.db_out_file)


    def create_gffdb(self):
        """Create a gffutils database from an ensembl gff3 file"""
        # logger.info("Generating DB from %s. (This may take 10-20+ minutes)", self.gff_in_file)

        gffutils.create_db(
            self.gff_in_file,
            self.db_out_file,
            force_gff=True,
            force=True,
            merge_strategy="create_unique"
        )
        # logger.info("DB generation complete [output_file=%s].", self.db_out_file)


    def test_db(self):
        """Test whether an existing database can be loaded and accessed"""
        # logger.info("Testing whether db can be accessed...")
        gff3_db = gffutils.FeatureDB(self.db_out_file)
        for feature in gff3_db.featuretypes():
            feature_counts = gff3_db.count_features_of_type(feature)
            # logger.info("Found %s features of type %s in the db.", feature_counts, feature)
            pass


    def run(self, overwrite: bool):
        """Run the script"""
        self.clean_overwrite(overwrite)
        self.create_gffdb()
        self.test_db()


def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description="Generate mapping stats for fusion detection")
    parser.add_argument(
        '--gff_input',
        dest='gff_in',
        help='ensembl gff3 annotation file',
        required=True
    )
    parser.add_argument(
        '--db_output',
        dest='db_out',
        help='gff3 derived annotation database from gffutils',
        required=True
    )
    parser.add_argument(
        '--overwrite',
        dest='overwrite',
        help='Overwrite existing database if present',
        default=False,
        action='store_true'
    )
    args = parser.parse_args()
    Gff2Db(args.gff_in, args.db_out).run(args.overwrite)


if __name__ == '__main__':
    main()
