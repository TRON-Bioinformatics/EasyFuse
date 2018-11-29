#!/usr/bin/env python

import os
import sys
import time
import subprocess

from argparse import ArgumentParser


def get_jobs_by_name(name, system="slurm"):
    '''This function returns a list of job ids for a given job name.'''
    if system == "slurm":
        return get_jobs_by_name_slurm(name)
    elif system == "pbs":
        return get_jobs_by_name_pbs(name)

def get_jobs_with_name_pbs(name):
    jobs = []
    output = ""
    try:
        output = subprocess.check_output(["qstat", "-f"])
    except:
        output = subprocess.Popen(["qstat", "-f"],stdout=subprocess.PIPE).communicate()[0]

    job_id = ""
    job_name = ""
    for x in output.rstrip().split("\n"):
        if x.lstrip().startswith("Job Id"):
            job_id = x.lstrip().rsplit(" ",1)[1]
        if x.lstrip().startswith("Job_Name"):
            job_name = x.lstrip().split(" = ")[1]

            if name == job_name:
                jobs.append(job_id)
            job_id = ""
            job_name = ""
    return jobs

def get_jobs_by_name_slurm(name):
    output = ""
    try:
        output = subprocess.check_output(["squeue", "-o %.18i %.9P %.60j %.8u %.2t %.10M %.6D %R"])
    except:
        output = subprocess.Popen(["squeue", "-o %.18i %.9P %.60j %.8u %.2t %.10M %.6D %R"], stdout=subprocess.PIPE).communicate()[0]
    jobs = []
    for x in output.rstrip().split("\n")[1:]:
        elements = x.lstrip().split()

        if name in elements[2]:
            jobs.append(elements[0])
    return jobs


def submit(job_name, cmd, cores, mem_usage, output_results_folder, dependencies, partitions, sched="slurm"):
    if sched == "slurm":
        submit_slurm(job_name, cmd, cores, mem_usage, output_results_folder, dependencies, partitions)
    elif sched == "pbs":
        submit_pbs(job_name, cmd, cores, mem_usage, output_results_folder, dependencies)
    else:
        submit_nonqueue(cmd, output_results_folder)
    
def submit_nonqueue(cmd, output_results_folder):
    print(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_results_folder, shell=True)
    (stdoutdata, stderrdata) = p.communicate()
    r = p.returncode
    if r != 0:
        print(stderrdata)
        sys.exit(1)

def submit_pbs(job_name, cmd, cores, mem_usage, output_results_folder, dependencies):
    p = subprocess.Popen('qsub', stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    (output, input) = (p.stdout, p.stdin)

    processors = "nodes=1:ppn={},mem={}gb,vmem={}gb".format(cores, mem_usage, mem_usage)
    error_file = os.path.join(output_results_folder, "error.log")
    output_file = os.path.join(output_results_folder, "output.log")

    job_string = """#!/bin/bash
#PBS -N {}
#PBS -l walltime=999:00:00
#PBS -l {}
#PBS -W depend=afterok:{}
#PBS -d {}
#PBS -e {}
#PBS -o {}

{}""".format(job_name, walltime, processors, ":".join(dependencies), output_results_folder, error_file, output_file, cmd)

    # Send job_string to qsub
    input.write(job_string)
    input.close()

    print(job_string)
        
    time.sleep(1)


def submit_slurm(job_name, cmd, cores, mem_usage, output_results_folder, dependencies, partitions):
    '''This function submits a predefined job with specific SBATCH parameters to the Slurm workload manager system.'''
    p = subprocess.Popen('sbatch', stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    (output, input) = (p.stdout, p.stdin)

    depend = ""
    if len(dependencies) != 0:
        depend = "#SBATCH --dependency=afterok:{}".format(":".join(dependencies))
    else:
        depend = ""
    error_file = os.path.join(output_results_folder, "error.log")
    output_file = os.path.join(output_results_folder, "output.log")

    job_string = """#!/bin/bash
#SBATCH -J {}
#SBATCH -p {}
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --cpus-per-task={}
#SBATCH --mem={}
#SBATCH --time=30-00:00:00
{}
#SBATCH --workdir={}
#SBATCH --error={}
#SBATCH --output={}

srun {}""".format(job_name, partitions, cores, int(mem_usage)*1000, depend, output_results_folder, error_file, output_file, cmd)

    input.write(job_string)
    input.close()

    # Generate sbatch script for your job in working dir
    outf = open(os.path.join(output_results_folder, "{}.sbatch".format(job_name)), "w")
    print(job_string)
    outf.write(job_string)
    outf.close()

def main():
    parser = ArgumentParser(description='Handle data streams')
    parser.add_argument('-j', '--job-name', dest='job_name', help='Specify the job name.', required=True )
    parser.add_argument('-s', '--script', dest='script', help='Specify the script to process.', required=True )
    parser.add_argument('-c', '--cores', dest='cores', type=int, help='Specify the number of cores to run the script with.', required=True )
    parser.add_argument('-m', '--memory', dest='memory', type=int, help='Specify the amount of memory to run your script with.', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Select the output folder.', required=True)
    parser.add_argument('-p', '--partitions', dest='partitions', help='Select the slurm partitions for this job.', default='allNodes')
    args = parser.parse_args()

    job_name = args.job_name
    cmd = args.script
    cores = args.cores
    mem_usage = args.memory
    output_results_folder = args.output
    dependencies = ""
    submit(job_name, cmd, cores, mem_usage, output_results_folder, dependencies, args.partitions)

if __name__ == '__main__':
    main()
