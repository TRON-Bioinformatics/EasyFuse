# load libraries ---------------------------------------------------------------

library(optparse)
library(tidyverse)
library(randomForest)

# Load functions ---------------------------------------------------------------
source("R/flat_table.R")
source("R/prediction_functions.R")

# Parse commandline arguments --------------------------------------------------
options(stringsAsFactors = FALSE)

argument_list <- list(
	make_option(c("-i", "--detected_fusions"), default="", help="Input list of detected fusions"), 
	make_option(c("-c", "--context_seq"), default="", help="Annotated context sequence list of detected fusions"), 
	make_option(c("-q", "--quantification"), default="", help="Quantification values for all fusion genes"), 
	make_option(c("-qc", "--qc_table"), default="", help="QC table file"), 
	make_option(c("-t", "--tool_state"), default="", help="tool state file"),
	make_option(c("-m", "--model_file"), default="", help="Path to .rds file containing the machine learning model"),
	make_option(c("-o", "--output"), default="", help="Final Output file for predicted fusion genes"),
)
opt <- parse_args(OptionParser(option_list=argument_list))

# check mandatory arguments
if(is.na(opt$detected_fusions) | opt$detected_fusions == "") {
	stop("Mandatory parameter \"detected fusions file\" missing. Aborting...")
}
if(is.na(opt$context_seq) | opt$context_seq == "") {
  stop("Mandatory parameter \"context sequence file\" missing. Aborting...")
}
if(is.na(opt$quantification) | opt$quantification == "") {
  stop("Mandatory parameter \"quantification file\" missing. Aborting...")
}
if(is.na(opt$qc_table) | opt$qc_table == "") {
  stop("Mandatory parameter \"fast qc file\" missing. Aborting...")
}
if(is.na(opt$tool_state) | opt$tool_state == "") {
  stop("Mandatory parameter \"tool state file\" missing. Aborting...")
}
if(is.na(opt$model_file) | opt$model_file == "") {
	print("Mandatory parameter \"model file\" missing. Aborting...")
}
if(is.na(opt$output) | opt$output == "") {
  stop("Mandatory parameter \"output file\" missing. Aborting...")
}

# #parameter in opt structure 
# opt$detected_fusions

# Read fusions -----------------------------------------------------------------

# read fusion with annotation and requantification to single table
fusion_flat <- read_fusion(opt$detected_fusions, 
                           opt$context_seq, 
                           opt$quantification, 
                           opt$qc_table_file, 
                           opt$tool_state)

# Prepare data for modeling ----------------------------------------------------
flat_fusion <- flat_fusion %>% 
  
  # replace NA's in read counts with 0
  replace_NA_in_reads() %>% 
  
  # convert categorial variables to factors
  prepare_for_modeling()


# Add predictions --------------------------------------------------------------

# read model from input file
model <- readr::read_rds(opt$model_file)

# add predictions
flat_fusions <- add_prediction(flat_fusions, model, 0.5)


# Write fusions to output file -------------------------------------------------
write_tsv(flat_fusions, opt$output)

