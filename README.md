# EasyFuse 

EasyFuse is a pipline for fusion gene detection from RNA-seq data.\
The pipeline is mainly written in python with some scripts being R.


## Installation

Fusion breakpoint prediction itself is currently not implemented in EasyFuse and the pipeline therefore depends on external fusion prediction tools.\
For simplicity we provide in the following an installation instruction for EasyFuse together with STAR-Fusion (ref) and Fusioncatcher (ref).

### Dependencies

 - Python (2.7.15)
 - Python modules:
    - pandas (0.24.0)
    - matplotlib (2.2.3)
    - seaborn (0.9.0)
    - pysam (0.15.2)
    - star (2.6.1b)
    - crossmap (0.2.7) (optional if liftover shall be included)
    - star-fusion (1.5.0)
    - biopython (1.73)
    - xlrd (1.0.0)
    - openpyxl (2.5.0a2)
    - bowtie2 (2.3.4.3)
    - bx-python (0.8.2)

Install python modules (we strongly recommend installation via conda):

  ```
  /path/to/conda/bin/conda create -n easyfuse python=2.7.15
  source /path/to/conda/bin/activate easyfuse
  conda install -c conda-forge pandas=0.24.0 matplotlib=2.2.3 seaborn=0.9.0 biopython=1.73 xlrd=1.0.0 openpyxl=2.5.0a2
  conda install -c bioconda pysam=0.15.2 star=2.6.1b star-fusion=1.5.0 bowtie2=2.3.4.3 bx-python=0.8.2 crossmap=0.2.7
  ```


 - R (>= 3.5.1)
 - R packages: 
    - optparse
    - tidyverse
    - randomForest
    - Biostrings
    - GenomicRanges
    - BSgenome    
    - bindrcpp

  Install packages within R by
  
  ```
  install.packages(c("optparse", "tidyverse", "randomForest", "Biostrings","BiocManager","BSgenome","optparse"))
  BiocManager::install("GenomicRanges") #bioconductor package
  ```
  

## Usage


### Start Fusion Prediction Pipeline

To start the fusion prediction pipeline on a specific sample the following python script has to
be executed with the given input parameters as command-line arguments.

```
processing.py \
    -i <test_sample_folder> \
    -o <working_dir> \
    -u <account_name> \
    -c <config_file> \
    -p <queueing_partition> \
    --tool_support 1 \
```

