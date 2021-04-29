# EasyFuse 

EasyFuse is a pipline for fusion gene detection from RNA-seq data. More detailed documentation is available in the [EasyFuse Wiki](https://github.com/TRON-Bioinformatics/EasyFuse/wiki).

A manuscript describing the method and performance evaluations is submitted for peer-review and publication.

## Installation

Fusion breakpoint prediction itself is currently not implemented in EasyFuse and the pipeline therefore depends on external fusion prediction tools.\
Prediction tools that have been implemented and tested within EasyFuse are listed under Tools. EasyFuse requires [STAR](https://github.com/alexdobin/STAR) for alignments. Additional alignment tools might be required depending on the external fusion prediction tools.
For simplicity we provide in the following an installation instruction for EasyFuse together with [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) and [Fusioncatcher](https://github.com/ndaniel/fusioncatcher). Detailed installation instructions can be found in the [EasyFuse Wiki](https://github.com/TRON-Bioinformatics/EasyFuse/wiki)

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


### Configuration

Before executing the pipeline some configuration files need to be adopted:

 - rename `build_env.sh.smaple` into `build_env.sh` and configure content. 
 - rename `config.py.smaple` into `config.py` and configure content.
 - rename `blacklist.txt.sample` into `blacklist.txt`.

## Usage

### Start Fusion Prediction Pipeline

To start the fusion prediction pipeline on a specific sample the following python script has to
be executed with the given input parameters as command-line arguments.

```
processing.py \
    -i <sample_folder> \
    -o <working_dir> 
```

### Example call on test dataset:

```
python processing.py -i test_case/SRR1659960_05pc_* -o test_easyfuse_1.3.4/
```
## Output

The output of EasyFuse is described in the wiki page [EasyFuse Output](https://github.com/TRON-Bioinformatics/EasyFuse/wiki/EasyFuse-Output).
