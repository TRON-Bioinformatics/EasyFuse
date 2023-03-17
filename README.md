# EasyFuse 

[![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/tron-bioinformatics/EasyFuse)](https://github.com/TRON-Bioinformatics/EasyFuse/releases)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/tronbioinformatics/easyfuse?label=docker)](https://hub.docker.com/r/tronbioinformatics/easyfuse)
[![License](https://img.shields.io/badge/license-GPLv3-green)](https://opensource.org/licenses/GPL-3.0)

EasyFuse is a pipeline to detect fusion transcripts from RNA-seq data with high accuracy.
EasyFuse uses five fusion gene detection tools, [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki), [InFusion](https://bitbucket.org/kokonech/infusion/src/master/), [MapSplice2](http://www.netlab.uky.edu/p/bioinfo/MapSplice2), [Fusioncatcher](https://github.com/ndaniel/fusioncatcher), and [SoapFuse](https://sourceforge.net/p/soapfuse/wiki/Home/) along with a powerful read filtering strategy, stringent re-quantification of supporting reads and machine learning for highly accurate predictions.

<p align="center"><img src="img/easyfuse_workflow.png"></p>

 - Publication: [Weber D, Ibn-Salem J, Sorn P, et al. Nat Biotechnol. 2022](https://doi.org/10.1038/s41587-022-01247-9)

We recommend using EasyFuse with Docker or Singularity.


## Installation

Fusion breakpoint prediction itself is currently not implemented in EasyFuse and the pipeline therefore depends on external fusion prediction tools.\
Prediction tools that have been implemented and tested within EasyFuse are listed under Tools. EasyFuse requires star for alignment, additional alignment
tools might be required depending on the external fusion prediction tools.
For simplicity we provide in the following an installation instruction for EasyFuse together with STAR-Fusion (ref) and Fusioncatcher (ref).

### Tools

 - star-fusion (1.5.0)
 - fusioncatcher (1.00)
 - mapsplice2 (2.2.1)
 - infusion (0.8)
 - SOAPfuse (1.2.7)
 - pizzly (0.37.3)
 - bowtie2 (2.3.4.3)
 - kallisto (0.44.0)
 - skewer (0.2.2)

### Dependencies

 - gzip (>=1.6) (if mapsplice shall be used)
 - samtools (1.9)
 - star (2.6.1d) 
 - Python (2.7.15)
 - Python modules:
    - pandas (0.24.0)
    - matplotlib (2.2.3)
    - seaborn (0.9.0)
    - pysam (0.15.2)
    - crossmap (0.2.7) (optional if liftover shall be included)
    - biopython (1.73)
    - xlrd (1.0.0)
    - openpyxl (2.5.0a2)
    - bx-python (0.8.2)

Install python modules (we strongly recommend installation via conda):

  ```
  /path/to/conda/bin/conda create -n easyfuse python=2.7.15
  source /path/to/conda/bin/activate easyfuse
  conda install -c conda-forge pandas=0.24.0 matplotlib=2.2.3 seaborn=0.9.0 biopython=1.73 xlrd=1.0.0 openpyxl=2.5.0a2
  conda install -c bioconda pysam=0.15.2 star=2.6.1b star-fusion=1.5.0 bowtie2=2.3.4.3 bx-python=0.8.2 crossmap=0.2.7
  ```

 - R (>= 3.6.0)
 - R packages: 
    - optparse (1.6.4)
    - tidyverse (1.3.0)
    - randomForest (4.6-14)

  Install packages within R by
  
  ```
  install.packages(c("optparse", "tidyverse", "randomForest"))
  ```
  
## Usage


### Start Fusion Prediction Pipeline

To start the fusion prediction pipeline on a specific sample the following python script has to
be executed with the given input parameters as command-line arguments.

```
processing.py \
    -i <test_sample_folder> \
    -o <working_dir> \
    --tool_support 1 \
```

### Example call (Test Case):

Before executing the example command

 - rename `build_env.sh.smaple` into `build_env.sh` and configure content. 
 - rename `config.py.smaple` into `config.py` and configure content.
 - rename `blacklist.txt.sample` into `blacklist.txt`.

```
python processing.py -i test_case/SRR1659960_05pc_* -o test_easyfuse_1.3.1/
```
## NextFlow

### Background
NextFlow is a workflow engine that allows to define and run workflows with multiple steps using a Domain Specific Language (DSL). It provides an expressive syntax that allows to parallelize tasks, scale out by parallelization across multiple machines seamlessly, define computational requirements with a high granularity, audit workflow performance, virtualise hardware and software requirements and improve reproducibility across different users and environments. Nextflow integrates with several queue batch processors and in particular with slurm. 
Installation and setup
NextFlow runs as a process for each workflow being executed, there is no need to setup any service.
Download the binary as follows:

```
curl -s https://get.nextflow.io | bash
```

You will need java >= 8. The default java in tronland is OpenJDK Java 8, but OpenJDK sometimes cause difficult to spot issues, I would recommend using the Oracle Java 11 that can be loaded as follows:
```
module load java/11.0.1
```
There is a common configuration file where you can define among other things the queue batch processor to use or the default amount of memory for each task.The configuration file has to be located at $HOME/.nextflow/config
To run nextflow in Tronland it is enough with the following configuration:
```
process.executor = 'slurm'

#this is optional

process.memory = '4G'
```
You can also use different profiles to switch between settings more easily. Below the standard profile is used by default, but if indicated with "-profile test" it will use the test profile and not run in slurm.
```
profiles {
    standard {
        process{
            executor = 'slurm'
            memory = '4G'
            clusterOptions = '-x tronc6,tronc7'
            errorStrategy = "finish"
        }
    }
    test {
        process.executor = 'local'
    }
}
```
It may be important not to use all nodes in the cluster to avoid saturating the cluster with something like "clusterOptions = '-x tronc6,tronc7'".

Also, when running multiple lengthy jobs in parallel and the first fails you may want that the rest to finish in order to cache their results; this can be achieved using "errorStrategy=finish".

It may be handy to setup nextflow so intermediate files are always deleted after execution. The ability to cache certain steps in a pipeline will be lost of course and you have to ensure that you use the directive "publishDir" with mode "copy" or "move" so your final outputs are not removed. Add this to your config file if you understand the implications.
```
cleanup = true
```
For further details on configuration see https://www.nextflow.io/docs/latest/config.html#configuration-file

## Sample workflow
The following is a sample process that receives a BAM file and a name prefix as input and outputs the deduplicated BAM with its index and a file of duplicate reads.
```
process picardDeduplication {
    module 'java/11.0.1'
    memory '16g'
    cpus 1
    publishDir "${output_folder}", mode: "copy"
    tag "${name}"
 
    input:
      file bam
      val name
 
    output:
      file "${name}.dedup.bam" into dedup_bam
      file "${name}.dedup.bai" into dedup_bai
      file "${name}.dedup.bam.bai" into dedup_bam_bai
      file "${name}.duplicates.txt" into duplicates
 
    """
    java -Djava.io.tmpdir=${scratch_folder} -Xmx${task.memory} -jar ${picard_jar} MarkDuplicates I=${bam} M=${name}.duplicates.txt O=${name}.dedup.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true CREATE_INDEX=true
    cp '${name}.dedup.bai' '${name}.dedup.bam.bai'
    """
}
```
NOTE: the publishDir defines where your outputs will be linked from, copied or moved  as defined. For large jobs it is important to define with care to avoid the later need to pick results from the miriad of working subfolders.

NOTE 2: the tag is relevant to define to identify individual jobs. For instance, if one given sample takes triple the time that all other samples for the same job you may want to identify which sample caused that. You can set as many tags as you wish.



This process has to be used in a workflow, as follows:
```
#!/usr/bin/env nextflow
 
/**
This workflow is a simple one step workflow that runs Picard's deduplication
*/
 
picard_jar = "/code/picard/2.18.17/picard.jar"
scratch_folder = "/scratch/info/projects/my_dummy_project"
 
params.name = false
params.bam = false
 
process picardDeduplication {
 
    [...]
 
}
```
We can include multiple processes in a workflow and if they do not depend on each, as in the example below, they will be parallelized:

```
#!/usr/bin/env nextflow
 
/**
This workflow is a simple one step workflow that runs Picard's deduplication
*/
 
picard_jar = "/code/picard/2.18.17/picard.jar"
scratch_folder = "/scratch/info/projects/my_dummy_project"
 
params.name = false
params.bam = false
 
process picardDeduplication {
 
    [...]
 
}
 
process snvVariantCalling {
    input:
        file "${name}.dedup.bam" into dedup_bam
        file "${name}.dedup.bai" into dedup_bai
 
    [...]
}
 
process indelVariantCalling {
    input:
        file "${name}.dedup.bam" into dedup_bam
        file "${name}.dedup.bai" into dedup_bai
 
    [...]
}
```

## Run a workflow
To run a workflow we just need to call the binary nextflow with the workflow file and the required parameters, such as:

```
nextflow my_picard_dedup_workflow.nf --bam my_input.bam --name my_deduplicated
```
This will trigger the orchestrator nextflow process which will send jobs to slurm for execution. The nextflow process needs to be running for the whole duration of the workflow execution. For this reason it is recommended to run the nextflow process inside a screen, a tmux or run it as job in the queue. The computational requirements for the nextflow process are minimal, run it as follows:

```
sbatch -N 1 --mem=1G -n 1 --job-name NF run_sample.sh
```

where in run_sample.sh you have:
```
#!/bin/bash
 
module load java/11.0.1
nextflow my_picard_dedup_workflow.nf --bam my_input.bam --name my_deduplicated
```
NOTE: the folder where you run nextflow will be where intermediate files and logs will be stored and depending on your workflow read and write operations will happen in this folder. It is recommended to use folders under /flash to profit from SSD read/write speed. This can be combined with a "publishDir" value in /scratch/info so you make sure that your results are stored in a more persistent file system.

WARNING: any module loaded in the session where nextflow is triggered from will be inherited by all processes potentially causing conflicts, ensure to run "module purge" before executing any nextflow workflow

Rexecuting workflows
Did your workflow failed in the third step after 10 hours of execution? No worries, just rerun it using the flag "-resume". Make sure that you don't mess with the intermediate files stored by nextflow, any modification will avoid you from rexecuting.

See the details here https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html

