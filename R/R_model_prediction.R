#!/usr/bin/env Rscript

# load libraries ---------------------------------------------------------------

library(optparse)
require(dplyr)
require(readr)
require(magrittr)
library(randomForest)

# Load functions ---------------------------------------------------------------

# Parse commandline arguments --------------------------------------------------
options(stringsAsFactors = FALSE)

argument_list <- list(
	make_option(c("-i", "--fusion_summary"), default="", help="Input list of detected fusions"), 
	make_option(c("-m", "--model_file"), default="", help="Path to .rds file containing the machine learning model"),
	make_option(c("-t", "--prediction_threshold"), default=0.75, help="Prediction score threshold above which a fusion is called positive"),
	make_option(c("-o", "--output"), default="", help="Final Output file for predicted fusion genes")
)
opt <- parse_args(OptionParser(option_list=argument_list))

# check mandatory arguments
if(is.na(opt$fusion_summary) | opt$fusion_summary == "") {
	stop("Mandatory parameter \"Input summary fusions file\" missing. Aborting...")
}
if(is.na(opt$model_file) | opt$model_file == "") {
	print("Mandatory parameter \"model file\" missing. Aborting...")
}
if(is.na(opt$output) | opt$output == "") {
  stop("Mandatory parameter \"output file\" missing. Aborting...")
}


# Read fusions -----------------------------------------------------------------

# read fusion with annotation and requantification to single table
fusion_data <- read_delim(opt$fusion_summary, delim = ";") %>%
  # convert character vectors into factors
  mutate(
    type = factor(type, levels = c("cis_near", "cis_inv", "cis_far", 
                                   "cis_trans", "trans", "trans_inv")),
    exon_boundary = factor(exon_boundary, levels = c("none", "3prime", 
                                                     "5prime", "both")),
    frame = factor(frame, levels = c("no_frame", "in_frame", "out_frame", 
                                     "neo_frame"))
  )

#' Function to add prediction columns according to imput model
#'
#' @param flat_data A data.frame with fusions.
#' @param model model object of type randomForest (See: randomForest package)
#' @param threshold single numeric as threshold. Events with prediction score >=
#'   threshold will be labled "positve" in prediction_class column and
#'   "negative" otherwise.
#'
#' @return a data.frame like input but with two additional columns named
#'   "prediction_prob" and "prediction_class".
#'   
add_prediction <- function(flat_data, model, threshold){
  
  # test if all features from the model are contained in the input data
  model_features <- rownames(model$importanceSD)
  model_features_contained <- all(model_features %in% names(flat_data))
  
  # if not all featrues appear in data, add only NAs
  if (! model_features_contained) {
    flat_data <- flat_data %>% 
      mutate(
        prediction_prob = NA,
        prediction_class = NA
      )
    
  # if all featrues appear in data add predictions
  } else {
    
    # get prediction probabilites
    pred_matrix <- predict(model, newdata = flat_data, type = "prob")
    pred_prob <- pred_matrix[, "positive"]
    
    # add columns
    flat_data <- flat_data %>% 
      mutate(
        prediction_prob = pred_prob,
        prediction_class = ifelse(pred_prob >= threshold, "positive", "negative")
      )
  }
  return(flat_data)
}


# read model from input file
model <- readr::read_rds(opt$model_file)

# add predictions
fusion_data <- add_prediction(fusion_data, model, opt$prediction_threshold)

# Write fusions to output file -------------------------------------------------
#consitently use dec="."
write.table(fusion_data,  opt$output, sep=";", dec=".", row.names = FALSE)

