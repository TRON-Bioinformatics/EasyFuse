#!/usr/bin/env python3

"""
Input/Output method collection
- Create folder and grant permission
- Read fastq files from directories
- get basename from full path

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""

import os
import re
import stat



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
    return fastqs

def pair_fastq_files(fastqs, logger=None):
    left = []
    right = []
    sample_id = []
    # iterate over the sorted list of file names and check for left/right pair file
    print("\nGoing to process the following read files...")
    for i, fq_file in enumerate(sorted(fastqs)):
        print("Fastq file {0}: {1}".format(i, fq_file))
        try:
            # Search for 1 or 2 between "_R|_" and "_|.f" in the basename of the file
            matches = re.search('(.*)(_R?[1-2])(_|[.])', os.path.basename(fq_file))

            sample_id = matches.group(1)
            forrev = matches.group(2)

            if forrev == "_1" or forrev == "_R1":
                left.append((sample_id, fq_file))
            elif forrev == "_2" or forrev == "_R2":
                right.append((sample_id, fq_file))


        except AttributeError:
            forrev = '-1'

    # Check whether file names match between paired files
    sample_list = []
    #if right:
    for sample_id_a, fq_file_a in left:
        found = False
        for sample_id_b, fq_file_b in right:
            if sample_id_a == sample_id_b:
                #print(sample_id_a, fq_file_a, fq_file_b)
                found = True
                sample_list.append((sample_id_a, fq_file_a, fq_file_b))
                break
            
        if not found:
            if logger:
                logger.error('Could not find a matching FASTQ for {}!'.format(sample_id_a))
            print('Could not find a matching FASTQ for {}!'.format(sample_id_a))
            sample_list.append((sample_id_a, fq_file_a, ""))

    return sample_list
