# EasyFuse 

EsayFuse is a pipline for fusion gene detection from RNA-seq data.


## Installation

### Dependencies

 - R (>= 3.5.1)
 - R packages: 
 
    - optparse
    - tidyverse
    - randomForest

  Install packages within R by
  
  ```
  install.packages(c("optparse", "tidyverse", "randomForest"))
  ```
  

## Usage

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

