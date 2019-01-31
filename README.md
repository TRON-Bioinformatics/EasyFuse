# EasyFuse 

EsayFuse is a pipline for fusion gene detection from RNA-seq data.


## Installation

### Dependencies

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
python processing.py \
    -i <test_sample_folder> \
    -o <working_dir> \
    -u <account_name> \
    -c <config_file> \
    -p <queueing_partition> \
    --tool_support 1 \
```


### Annotate fusion breakpoints

To annotate predicted breakpoints for fusion genes the following R script has to 
be executed with the given input files as command-line arguments.

```
Rscript R/GetFusionSequence.R \
	-i <detected_fusion_file> \
	-f <fasta_genome_dir> \
	-e <ensemble_csv_file> \
	-o <output_file> \
	--cis_near_distance \
	--genomic_seq_len \
	--context_seq_len \
```	

### Apply prediction model

To apply the prediction model to detected fusions the following R script has to 
be executed with the given input files as command-line arguments.

```
Rscript R/model_prediction.R \
	-i <detected_fusion_file> \
	-c <context_seq_file> \
	-q <requantification_file> \
	-qc <qc_table_file> \
	-t <tool_state_file> \
	-m <model_file> \
	-o <output_file> 
```

