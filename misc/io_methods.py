#!/usr/bin/env python

"""
Input/Output method collection
- Create folder and grant permission
- Read fastq files from directories
- get basename from full path

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""

from __future__ import print_function
import os
import sys
import stat
import re
import ntpath


def create_folder(path, logger=None):
    '''This function creates a specified folder, if not already existing and grants the respective permission.'''
    if not os.path.exists(path):
        os.makedirs(path)
        grant_permissions(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)
        if logger:
            logger.debug("Created folder {}".format(path))
        else:
            print("Created folder {}".format(path))
    else:
        if logger:
            logger.debug("Folder {} already exists".format(path))
        else:
            print("Folder {} already exists".format(path))


def grant_permissions(path, mode):
    '''This function grants specific permissions for a given folder or file.'''
    for _, dirs, files in os.walk(path, topdown=False):
        for dir in dirs:
            os.chmod(dir, mode)
        for file in files:
            os.chmod(file, mode)


def get_fastq_files(input_paths, logger=None):
    """Load fastq files, check file name consistency and return tuple of left/right pairs"""
    # load fastq files from the provide path to a list
    fastqs = []
    for path in input_paths:
        if os.path.isdir(path):
            files = os.listdir(path)
            for filename in files:
                file_path = os.path.join(path, filename)
                if os.path.isfile(file_path) and filename.endswith((".fastq.gz", ".fq.gz", ".fq", ".fastq")):
                    if not filename.endswith(".gz"):
                        logger.info("Warning: Fastq files are supposed to be gzipped! Skipping record {}!".format(file))
                        continue
                    fastqs.append(file_path)
                elif os.path.isdir(file_path):
                    fastqs_tmp = get_fastq_files([file_path])
                    fastqs.extend(fastqs_tmp)
        elif os.path.isfile(path) and path.endswith((".fastq.gz", ".fq.gz", ".fq", ".fastq")):
            fastqs.append(path)

    left = []
    right = []
    sample_id = []
    # iterate over the sorted list of file names and check for left/right pair file
    print("\nGoing to process the following read files...")
    for i, fq_file in enumerate(sorted(fastqs)):
        print("Fastq file {0}: {1}".format(i, fq_file))
        # urla: handled cases are: abc_R1_xyz.fastq.gz || abc_R1_xyz.fq.gz || abc_R1.fastq.gz || abc_1.fastq.gz || etc
        try:
            # Search for 1 or 2 between "_R|_" and "_|.f" in the basename of the file
            #forrev = re.search('(\_|\_R)([1-2])(\_|\.f)', path_leaf(fq_file)).group(2)
            # urla: replaced original (up) with the following <- worked in minimal test set, should be fine for everything
            forrev = re.search('_R?([1-2])(_|[.]f)', path_leaf(fq_file)).group(1)
#            print(forrev)
        except AttributeError:
            forrev = '-1'
                                
        if forrev == '1':
            left.append(fq_file)
        elif forrev == '2':
            right.append(fq_file)
        else:
            print('Warning: Ignoring \"{}\" as it doesn\'t seem to be a valid fastq file'.format(fq_file))
    # Check whether file names match between paired files
    if right:
        for i, _ in enumerate(left):
            # urla: ids for comparison are not generally applicable here, as they remove lane numbers as well
            #sample_id_l = re.search('(\w*)(\_|\_R)([1])(\_|\.f)', path_leaf(left[i])).group(1)
            #sample_id_r = re.search('(\w*)(\_|\_R)([2])(\_|\.f)', path_leaf(right[i])).group(1)
            # urla: replaced original (up) with the following <- worked in minimal test set, should be fine for everything
            sample_id_l = re.search('(.*)_R?1(_|[.]f)', path_leaf(left[i])).group(1)
            sample_id_r = re.search('(.*)_R?2(_|[.]f)', path_leaf(right[i])).group(1)
            sample_id.append(sample_id_l)
#            print(sample_id_l)
#            print(sample_id_r)
            if sample_id_l != sample_id_r:
                if logger:
                    logger.error('Error 99: Paired files names {0} and {1} do not match!'.format(sample_id_l, sample_id_r))
                print('Error 99: Paired files names {0} and {1} do not match!'.format(sample_id_l, sample_id_r))
                sys.exit(99)
    return (left, right, sample_id)


def path_leaf(path):
    """Return the basename of a file path"""
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
