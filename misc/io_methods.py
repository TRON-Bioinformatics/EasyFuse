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
import pwd
import grp
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
    for root, dirs, files in os.walk(path, topdown=False):
        for dir in dirs:
            os.chmod(dir, mode)
        for file in files:
            os.chmod(file, mode)


def get_fastq_files(input_path, logger=None):
    """Load fastq files, check file name consistency and return tuple of left/right pairs"""
    # load fastq files from the provide path to a list
    fastqs = []
    if os.path.isdir(input_path):
        files = os.listdir(input_path)
        for filename in files:
            file_path = os.path.join(input_path, filename)
            if os.path.isfile(file_path) and filename.endswith((".fastq.gz", ".fq.gz", ".fq", ".fastq")):
                if not filename.endswith(".gz"):
                    if logger:
                        logger.info("Warning: Fastq files are supposed to be gzipped! Skipping record {}!".format(file))
                    else:
                        print("Warning: Fastq files are supposed to be gzipped! Skipping record {}!".format(file))
                    continue
                fastqs.append(file_path)
            elif os.path.isdir(file_path):
                if logger:
                    logger.info("Not following detected nested dir \"{0}\" within fastq path input \"{1}\"".format(file_path, input_path))
                else:
                    print("Not following detected nested dir \"{0}\" within fastq path input \"{1}\"".format(file_path, input_path))
    elif os.path.isfile(input_path):
        if logger:
            logger.error("Error 99: File input is not supported, please provide a directory containing fastq files to process!")
        print("Error 99: File input is not supported, please provide a directory containing fastq files to process!")
        sys.exit(99)
                        
    left = []
    right = []
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
        for i, fq_file in enumerate(left):
            # urla: ids for comparison are not generally applicable here, as they remove lane numbers as well
            #sample_id_l = re.search('(\w*)(\_|\_R)([1])(\_|\.f)', path_leaf(left[i])).group(1)
            #sample_id_r = re.search('(\w*)(\_|\_R)([2])(\_|\.f)', path_leaf(right[i])).group(1)
            # urla: replaced original (up) with the following <- worked in minimal test set, should be fine for everything
            sample_id_l = re.search('(.*)_R?1(_|[.]f)', path_leaf(left[i])).group(1)
            sample_id_r = re.search('(.*)_R?2(_|[.]f)', path_leaf(right[i])).group(1)
#            print(sample_id_l)
#            print(sample_id_r)
            if sample_id_l != sample_id_r:
                if logger:
                    logger.error('Error 99: Paired files names {0} and {1} do not match!'.format(sample_id_l, sample_id_r))
                print('Error 99: Paired files names {0} and {1} do not match!'.format(sample_id_l, sample_id_r))
                sys.exit(99)
    return (left, right, sample_id_l)
                    

def get_icam_reads_data(root_path, sample_file_keys, logger):
    """Return tuple of lists of demultiplexing stats, and fastq dirs"""
    #root_path = "/mnt/bfx/IVAC2_0/BNT/iCaM2Scheduler/data/demux_out"
    # grep -A4 "Patient ID" /mnt/bfx/IVAC2_0/BNT/iCaM2Scheduler/data/demux_out/*/demux_stats.csv | cut -d";" -f-3 > demuxStats.txt
    sample_readcnt_dict = {}
    sample_fastq_dict = {}
    counter_list = [0, 0, 0, 0, 0] # 0:processed, 1:new, 2:unexpected format, 3: low read count, 4: fastq file error
    is_sample_ok = True
    for root, _, files in os.walk(root_path):
        # search "demux_stats.csv" to get patient id and directory/file names of fastq files
        for file_name in files:
            if file_name == "demux_stats.csv":
                is_sample_ok = True
                with open(os.path.join(root, file_name), "r") as dmx_file:
                    count_lines = 0
                    sample_rep1_id = ""
                    sample_rep2_id = ""
                    for line in dmx_file:
                        count_lines += 1
                        line_splitter = line.split(";")
                        # check for specific format, ignore everything else as it would introduce severe errors in downstream processing
                        # urla - todo: define rules for other formats (and why they are there)
                        if count_lines == 2 and (line_splitter[0] != "Patient ID" or line_splitter[1] != "Tumor RNA Reads Rep1" or line_splitter[2] != "Tumor RNA Reads Rep2"):
                            logger.debug("Unexpected input format, skipping file {} in order to avoid parsing errors".format(os.path.join(root, file_name)))
                            counter_list[2] += 1
                            break
                        # parse desired content only
                        if count_lines in [3, 5] and len(line_splitter[0]) > 1:
                            sample_rep1_id = "{}_Rep1_{}".format(line_splitter[0], line_splitter[1])
                            sample_rep2_id = "{}_Rep2_{}".format(line_splitter[0], line_splitter[2])
                        if count_lines in [4, 6] and len(line_splitter[0]) > 1:
                            if int(line_splitter[1]) < 75000000:
                                logger.debug("Found less than 75M reads, skipping sample {} due to likely errors".format(sample_rep1_id))
                                counter_list[3] += 1
                                is_sample_ok = False
                            else:
                                sample_readcnt_dict[sample_rep1_id] = int(line_splitter[1])
                                
                            if int(line_splitter[2]) < 75000000:
                                logger.debug("Found less than 75M reads, skipping sample {} due to likely errors".format(sample_rep2_id))
                                counter_list[3] += 1
                                is_sample_ok = False
                            else:
                                sample_readcnt_dict[sample_rep2_id] = int(line_splitter[2])
                # search for fastq files in the subdir of the sample folder
                # the name of the fastq dir is constructed from the number of identified reads
                for root2, _, files2 in os.walk(root):
                    for sample in sample_readcnt_dict:
                        if "Sample_{}".format(sample.split("_")[2]) in root2:
                            if len(files2) != 2:
                                logger.debug("Couldn't identify correct fastq files for sample {}. Skipping this one...".format(sample))
                                counter_list[4] += 1
                                is_sample_ok = False
                            else:
                                if "_R2_" in files2[0]:
                                    files2_tmp = files2[0]
                                    files2[0] = files2[1]
                                    files2[1] = files2_tmp
                                if is_sample_ok:
                                    sample_fastq_dict[sample] = (os.path.join(root, root2, files2[0]), os.path.join(root, root2, files2[1]))
                                if "_".join(sample.split("_")[:2]) in sample_file_keys:
                                    counter_list[0] += 1
                                else:
                                    counter_list[1] += 1
                            #print("read files? {}".format(files2))
    return (sample_readcnt_dict, sample_fastq_dict, counter_list)


# urla: this method is currently not used anywhere, but provide useful info in later releases
#       e.g. if the tumor content is rather low, it is more difficult to identify the fusion event in the wildtyp background;
#       vice versa, if the tumor content is rather high, a fusion event should be "easy" to find and have "stronger" support
def get_tumor_content_estimate(root_path, sample_id):
    """Return the tumor content derived from the log file of MyMut"""
    # grep "INFO est. Tumor content" /mnt/bfx/IVAC2_0/BNT/iCaM2Scheduler/data/iCaM2Bot_out/*/*.targets.confirmed.txt > tumorContent.txt
    root_path = "/mnt/bfx/IVAC2_0/BNT/iCaM2Scheduler/data/iCaM2Bot_out"
    sample_id = "A33438"
    for root, _, files in os.walk(root_path):
        for file_name in files:
            if "{}.targets.confirmed.txt".format(sample_id) in file_name:
                with open(os.path.join(root, file_name), "r") as tsc_file:
                    for line in tsc_file:
                        if "INFO est. Tumor content" in line:
                            print("file: {}; line: {}".format(file_name, line))
                            
    tumor_content = ""
    return tumor_content


def path_leaf(path):
    """Return the basename of a file path"""
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
