#!/usr/bin/env python

"""
Print done to the log file
This script is called after summarize_data finished successfully
It is a separate script a allow queuing compatibility

@author: BNT (URLA)
@version: 20190613
"""
from __future__ import print_function
import argparse
from logger import Logger

# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class AllDone(object):
    """Run, monitor and schedule fastq processing for fusion gene prediction"""
    def __init__(self, logger_path):
        """Parameter initiation and work folder creation. Start of progress logging."""
        self.logger = Logger(logger_path)
        self.logger.info("done")

def main():
    """Parse command line arguments and start script"""
    parser = argparse.ArgumentParser(description='Write final done if everything is done :)')
    # required arguments
    parser.add_argument('--logger', dest='logger', help='Path to the log file', required=True)
    args = parser.parse_args()
    AllDone(args.logger)

if __name__ == '__main__':
    main()
