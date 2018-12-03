"""
Simple method collection for:
1) Job submission to slurm workload manager (Creates an sbatch file and starts it)
2) Program execution outside the scheduler with "subprocess"
3) Slurm queue querying

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""
from __future__ import print_function
import os
import sys
import subprocess

# pylint: disable=line-too-long
#         yes they are partially, but I do not consider this to be relevant here
class Queue(object):
    """Prepare and run a job either via slurm or not"""
    def __init__(self):
        pass

    @staticmethod
    def submit_nonqueue(cmd):
        """Run a command directly on the console"""
        proc = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=False)
        stderrdata = proc.communicate()
        proc_return = proc.returncode
        if proc_return != 0:
            print("Error: Command \"{}\" returned non-zero exit status".format(cmd))
            print(stderrdata)
            sys.exit(100)

    @staticmethod
    def submit_nonqueue_and_get_stdout(cmd):
        """Run a command directly on the console and return the stdout"""
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)
        stdoutdata = proc.communicate()
        proc_return = proc.returncode
        if proc_return != 0:
            print("Error: Command \"{}\" returned non-zero exit status".format(cmd))
            sys.exit(100)
        return stdoutdata

    # pylint: disable=too-many-arguments
    #         all arguments are required for proper resource management
    @staticmethod
    def submit_slurm(job_name, cmd, cores, mem_usage, output_results_folder, dependencies, partitions, userid):
        """Submit a predefined job with specific SBATCH parameters to the Slurm workload manager system"""
        # Customize your options here
        depend = "\n"
        if len(dependencies) != 0: # dependencies is an array, not a str => pylint: disable=len-as-condition
            depend = "#SBATCH --dependency=afterok:{}\n".format(":".join(dependencies))
        error_file = os.path.join(output_results_folder, "error.log")
        output_file = os.path.join(output_results_folder, "output.log")

        # write the slurm sbatch script
        with open(os.path.join(output_results_folder, job_name + ".sbatch"), "w") as slurm_script:
            slurm_script.writelines([
                "#!/bin/bash\n",
                "#SBATCH -J {}\n".format(job_name),
                "#SBATCH -p {}\n".format(partitions),
                "#SBATCH --account {}\n".format(userid),
                "#SBATCH --kill-on-invalid-dep=yes\n",
                "#SBATCH --cpus-per-task={}\n".format(cores),
                "#SBATCH --mem={}\n".format(int(mem_usage)*1000),
                "#SBATCH --time=06:00:00\n",
                depend,
                "#SBATCH --workdir={}\n".format(output_results_folder),
                "#SBATCH --error={}\n".format(error_file),
                "#SBATCH --output={}\n".format(output_file),
                "\n",
                "##SBATCH --mail-type=ALL\n",
                "##SBATCH --mail-user=urs.lahrmann@biontech.de\n",
                "\n",
                "set -eo pipefail -o nounset\n",
                "srun {}\n".format(cmd)
                ])
        # run script
        try:
            output = subprocess.check_output(["sbatch", os.path.join(output_results_folder, job_name + ".sbatch")])
            print("{} ({})\n".format(output.rstrip(), job_name))
        except subprocess.CalledProcessError as call_error:
            print(call_error.output)
            print("Could not start slurm script {}".format(os.path.join(output_results_folder, job_name + ".sbatch")))
            sys.exit(100)

    @staticmethod
    def get_jobs_with_name_slurm(name):
        """Check if slurm job is already running (by name) and return its jobid if it does"""
        output = ""
        try:
            # urla: original was default output, simplified it to only JOBID and NAME, where NAME is being checked and JOBID returned
            output = subprocess.check_output(["squeue", "-o %i %j"])
        except subprocess.CalledProcessError as call_error:
            print(call_error.output)
            output = subprocess.Popen(["squeue", "-o %i %j"], stdout=subprocess.PIPE).communicate()[0]
        jobs = []
        # Check for every line of the squeue output if it contains "name"
        for job_list in output.rstrip().split("\n")[1:]:
            elements = job_list.lstrip().split()
            if name in elements[1]:
                jobs.append(elements[0])
        return jobs
