#!/usr/bin/env python3

"""
Simple method collection for:
1) Job submission to slurm workload manager (Creates an sbatch file and starts it)
2) Program execution outside the scheduler with "subprocess"
3) Slurm queue querying

@author: Tron (PASO), BNT (URLA)
@version: 20180118
"""

import os
import subprocess
import sys
from argparse import ArgumentParser

from logzero import logger


def get_jobs_by_name(name, system="slurm"):
    """This function returns a list of job ids for a given job name."""
    if system == "slurm":
        return get_jobs_by_name_slurm(name)
    elif system == "pbs":
        return get_jobs_by_name_pbs(name)
    else:
        return []


def get_jobs_by_name_pbs(name):
    jobs = []
    output = ""
    try:
        output = subprocess.check_output(["qstat", "-f"], universal_newlines=True)
    except:
        output = subprocess.Popen(
            ["qstat", "-f"], stdout=subprocess.PIPE
        ).communicate()[0]

    job_id = ""
    job_name = ""
    for x in output.rstrip().split("\n"):
        line = x.lstrip()
        if line.startswith("Job Id"):
            job_id = line.rsplit(" ", 1)[1]
        if line.startswith("Job_Name"):
            job_name = line.split(" = ")[1]

            if name == job_name:
                jobs.append(job_id)
            job_id = ""
            job_name = ""
    return jobs


def get_jobs_by_name_slurm(name):
    """Check if slurm job is already running (by name) and return its jobid if it does"""
    output = ""
    try:
        output = subprocess.check_output(
            ["squeue", "-o %i %j"], universal_newlines=True
        )
    except subprocess.CalledProcessError as call_error:
        logger.error(call_error.output)
        output = subprocess.Popen(
            ["squeue", "-o %i %j"], stdout=subprocess.PIPE
        ).communicate()[0]
    jobs = []
    for job_list in output.rstrip().split("\n")[1:]:
        elements = job_list.lstrip().split()

        if name in elements[1]:
            jobs.append(elements[0])
    return jobs


def submit(
    job_name,
    cmd,
    cores,
    mem_usage,
    output_results_folder,
    dependencies,
    partitions,
    userid,
    timelimit,
    mail,
    module_file,
    sched="slurm",
):
    if sched == "slurm":
        _submit_slurm(
            job_name,
            cmd,
            cores,
            mem_usage,
            output_results_folder,
            dependencies,
            partitions,
            userid,
            timelimit,
            mail,
            module_file,
        )
    else:
        _submit_nonqueue(job_name, cmd)


def _submit_nonqueue(job_name, cmd):
    #    if module_file:
    #        cmd = " && ".join(["source " + module_file, " ".join(cmd)]).split(" ")
    logger.info("Running {}".format(job_name))
    logger.info("CMD: {}".format(" ".join(cmd)))
    if ">" in cmd:
        index = cmd.index(">")
        output_file = cmd[index + 1]
        adj_cmd = cmd[:index]
        with open(output_file, "w") as fp:
            subprocess.run(adj_cmd, stdout=fp)
    else:
        p = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False
        )
        (_, stderrdata) = p.communicate()
        r = p.returncode
        if r != 0:
            logger.error(
                'Error: Command "{}" returned non-zero exit status'.format(cmd)
            )
            logger.error(stderrdata)
            sys.exit(1)


def _submit_slurm(
    job_name,
    cmd,
    cores,
    mem_usage,
    output_results_folder,
    dependencies,
    partitions,
    userid,
    timelimit,
    mail,
    module_file,
):
    """This function submits a predefined job with specific SBATCH parameters to the Slurm workload manager system."""
    # add dependencies if necessary
    depend = "\n"
    if dependencies:
        depend = "#SBATCH --dependency=afterok:{}\n".format(":".join(dependencies))
    # add mail delivery
    mail_type = "\n"
    mail_user = "\n"
    if mail:
        mail_type = "#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80\n"
        mail_user = "#SBATCH --mail-user={}\n".format(mail)
    # save stderr and stdout in different files
    error_file = os.path.join(output_results_folder, "error.log")
    output_file = os.path.join(output_results_folder, "output.log")

    # Generate sbatch script for your job in working dir
    slurm_script = os.path.join(output_results_folder, job_name + ".sbatch")
    with open(slurm_script, "w") as sbatch:
        sbatch.writelines(
            [
                "#!/bin/bash\n",
                "#SBATCH -J {}\n".format(job_name),
                "#SBATCH -p {}\n".format(partitions),
                "#SBATCH --account {}\n".format(userid),
                "#SBATCH --kill-on-invalid-dep=yes\n",
                "#SBATCH --cpus-per-task={}\n".format(cores),
                "#SBATCH --mem={}G\n".format(mem_usage),
                "#SBATCH --time={}\n".format(timelimit),
                depend,
                "#SBATCH -D {}\n".format(output_results_folder),
                "#SBATCH --error={}\n".format(error_file),
                "#SBATCH --output={}\n".format(output_file),
                mail_type,
                mail_user,
                "set -eo pipefail -o nounset\n",
                ". {}\n".format(module_file) if module_file else "\n",
                'srun echo "$(date) Slurm job started." > {0}.runtime && {1} && echo "$(date) Slurm job finished." >> {0}.runtime\n'.format(
                    job_name, cmd
                ),
            ]
        )
    # and run it
    try:
        output = subprocess.check_output(["sbatch", slurm_script])
        logger.info("{} ({})".format(output.rstrip(), job_name))
    except subprocess.CalledProcessError as call_error:
        logger.error(call_error.output)
        logger.error("Could not start slurm script {}".format(slurm_script))
        sys.exit(1)


def main():
    parser = ArgumentParser(description="Handle data streams")
    parser.add_argument(
        "-j", "--job-name", dest="job_name", help="Specify the job name.", required=True
    )
    parser.add_argument(
        "-s",
        "--script",
        dest="script",
        help="Specify the script to process.",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--cores",
        dest="cores",
        type=int,
        help="Specify the number of cores to run the script with.",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--memory",
        dest="memory",
        type=int,
        help="Specify the amount of memory to run your script with.",
        required=True,
    )
    parser.add_argument(
        "-o", "--output", dest="output", help="Select the output folder.", required=True
    )
    parser.add_argument(
        "-p",
        "--partitions",
        dest="partitions",
        help="Select the slurm partitions for this job.",
        default="allNodes",
    )
    parser.add_argument(
        "-u",
        "--userid",
        dest="userid",
        help="Select the slurm account to run the job on.",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--timelimit",
        dest="timelimit",
        help="Select the timelimit for the job.",
        default="30-00:00:0",
    )
    args = parser.parse_args()

    job_name = args.job_name
    cmd = args.script
    cores = args.cores
    mem_usage = args.memory
    output_results_folder = args.output
    dependencies = ""
    submit(
        job_name,
        cmd,
        cores,
        mem_usage,
        output_results_folder,
        dependencies,
        args.partitions,
        args.userid,
        args.timelimit,
        "",
        "",
    )
