#!/usr/bin/env python


import sys

from misc.samples import SamplesDB

def main():
    infile = sys.argv[1]

    db = SamplesDB(infile)

    db.print_db()


if __name__ == "__main__":
    main()
