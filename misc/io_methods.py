#!/usr/bin/env python

import os
import pwd
import grp
import stat


def create_folder(path):
    '''This function creates a specified folder, if not already existing and grants the respective permission.'''
    if not os.path.exists(path):
        print("Creating folder", path)
        os.makedirs(path)


def get_fastq_files(paths):
    '''This function returns a list of fastq files for the given list of paths.'''
    fastqs = []
    for path in paths:
        if os.path.isdir(path):
            files = os.listdir(path)
            for file in files:
                file_path = os.path.join(path,file)
                if os.path.isfile(file_path) and file.endswith((".fastq.gz","fastq")):
                    fastqs.append(file_path)
                elif os.path.isdir(file_path):
                    fastqs_tmp = IOMethods.get_fastq_files([file_path])
                    fastqs.extend(fastqs_tmp)
        elif os.path.isfile(path) and path.endswith((".fastq.gz","fastq")):
            fastqs.append(path)
    return fastqs


def grant_permissions(path, mode):
    '''This function grants specific permissions for a given folder or file.'''
    for root, dirs, files in os.walk(path, topdown=False):
        for dir in dirs:
            os.chmod(dir, mode)
        for file in files:
            os.chmod(file, mode)
