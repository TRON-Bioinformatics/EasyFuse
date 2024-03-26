#!/usr/bin/env Rscript

# load libraries ---------------------------------------------------------------
library(optparse, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(readr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(tidyselect, quietly = TRUE)
library(stringr, quietly = TRUE)
library(randomForest, quietly = TRUE)

# Parse commandline arguments --------------------------------------------------

argument_list <- list(
	make_option(c("-i", "--fusion_summary"), default="", 
	            help="Input list of detected fusions"), 
	make_option(c("-m", "--model_file"), default="", 
	            help="Path to .rds file containing the machine learning model"),
	make_option(c("-t", "--prediction_threshold"), default=0.5, 
	            help="Prediction score threshold above which a fusion is called 
	            positive"),
	make_option(c("-o", "--output"), default="", 
	            help="Final Output file for predicted fusion genes")
)
opt <- parse_args(OptionParser(option_list=argument_list))

# check mandatory arguments
if(is.na(opt$fusion_summary) | opt$fusion_summary == "") {
	stop("Mandatory parameter \"Input summary fusions file\" missing. 
	     Aborting...")
}
if(is.na(opt$model_file) | opt$model_file == "") {
	print("Mandatory parameter \"model file\" missing. Aborting...")
}
if(is.na(opt$output) | opt$output == "") {
  stop("Mandatory parameter \"output file\" missing. Aborting...")
}

input_file <- opt$fusion_summary
model_file <- opt$model_file
prediction_threshold <- opt$prediction_threshold
output_file <- opt$output

# Read fusion gene data --------------------------------------------------------

# read fusion with annotation and re-quantification to single table
fusion_data <- readr::read_delim(input_file, 
                                 delim = ";", 
                                 show_col_types = FALSE) %>%
  
  # convert character vectors into factors with fixed levels to match expected
  # types of the pre-trained model
  dplyr::mutate(
    
    type = factor(type, 
                  levels = c("cis_near", "cis_inv", "cis_far", "cis_trans", 
                             "trans", "trans_inv")),
    # fix NA in exon_boundary if present
    exon_boundary = tidyr::replace_na(exon_boundary, "no_match"),
    
    exon_boundary = factor(exon_boundary, 
                           levels = c("no_match", "3prime", "5prime", 
                                      "both")),
    frame = factor(frame, 
                   levels = c("no_frame", "in_frame", "out_frame", 
                              "neo_frame")),
  )

# read model from input file ---------------------------------------------------
model <- readr::read_rds(model_file)

# check if all features from the model are contained in the input data
model_features <- rownames(model$importance)
model_features_contained <- all(model_features %in% names(fusion_data))

# if not all features appear in data, add only NAs and raise a warning
if (! model_features_contained) {
  warning(
    str_c("Model features are not completely contained in input data. 
    Missing feature: ", 
    setdiff(model_features, names(fusion_data)))
  )
   fusion_data <- fusion_data %>%
    dplyr::mutate(
      prediction_prob = NA,
      prediction_class = NA
    )
  
# if all features appear in data add predictions
} else {
  
  # get prediction probabilities by appling the model to the input data
  pred_matrix <- predict(
    model, 
    newdata = dplyr::select(fusion_data, tidyselect::one_of(model_features)), 
    type = "prob")
  pred_prob <- pred_matrix[, "positive"]
  
  # add columns to input data
  fusion_data <- fusion_data %>%
    dplyr::mutate(
      prediction_prob = pred_prob,
      prediction_class = ifelse(
        pred_prob >= prediction_threshold, "positive", "negative")
    )
}

# Write fusions to output file -------------------------------------------------

# write only positive predicted fusion genes to output
readr::write_delim(
  dplyr::filter(fusion_data, prediction_class == "positive"), 
  output_file, 
  delim = ";")

# write additional output file with all candidate fusion genes
readr::write_delim(
  fusion_data, 
  paste0(stringr::str_replace(output_file, "\\.csv", ""), ".all.csv"), 
  delim = ";")
