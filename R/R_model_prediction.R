#!/usr/bin/env Rscript

# load libraries ---------------------------------------------------------------

require(optparse)
require(dplyr)
require(readr)
require(magrittr)
require(randomForest)

# Load functions ---------------------------------------------------------------

# Parse commandline arguments --------------------------------------------------
options(stringsAsFactors = FALSE)

argument_list <- list(
	make_option(c("-i", "--fusion_summary"), default="", help="Input list of detected fusions"), 
	make_option(c("-m", "--model_file"), default="", help="Path to .rds file containing the machine learning model"),
	make_option(c("-t", "--prediction_threshold"), default=0.5, help="Prediction score threshold above which a fusion is called positive"),
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

input_file <- opt$fusion_summary
model_file <- opt$model_file
prediction_threshold <- opt$prediction_threshold
output_file <- opt$output

# Read fusions -----------------------------------------------------------------

# read fusion with annotation and requantification to single table
fusion_data <- read_delim(input_file, delim = ";") %>%
  
  mutate(
    
    # convert character vectors into factors with fixed levels
    type = factor(type, 
                  levels = c("cis_near", "cis_inv", "cis_far", "cis_trans", 
                             "trans", "trans_inv")),
    # fix NA in exon_boundary if present
    exon_boundary = replace_na(exon_boundary, "no_match"),
    
    exon_boundary = factor(exon_boundary, 
                           levels = c("no_match", "3prime", "5prime", 
                                      "both")),
    frame = factor(frame, 
                   levels = c("no_frame", "in_frame", "out_frame", 
                              "neo_frame"))
  )

#' Function to add prediction columns according to imput model
#'
#' @param df A data.frame with fusions.
#' @param model model object of type randomForest (See: randomForest package)
#' @param threshold single numeric as threshold. Events with prediction score >=
#'   threshold will be labled "positve" in prediction_class column and
#'   "negative" otherwise.
#'
#' @return a data.frame like input but with two additional columns named
#'   "prediction_prob" and "prediction_class".
#'
add_prediction <- function(df, model, threshold){
  
  # test if all features from the model are contained in the input data
  model_features <- rownames(model$importance)
  model_features_contained <- all(model_features %in% names(df))
  
  # if not all featrues appear in data, add only NAs
  if (! model_features_contained) {
    warning("Model features are not completely contained in input data.")
    df <- df %>%
      mutate(
        prediction_prob = NA,
        prediction_class = NA
      )
    
    # if all featrues appear in data add predictions
  } else {
    
    # get prediction probabilites
    in_data <- select(df, all_of(model_features))
    pred_matrix <- predict(model, newdata = in_data, type = "prob")
    pred_prob <- pred_matrix[, "positive"]
    
    # add columns
    df <- df %>%
      mutate(
        prediction_prob = pred_prob,
        prediction_class = ifelse(
          pred_prob >= threshold, "positive", "negative")
      )
  }
  return(df)
}

# read model from input file
model <- readr::read_rds(model_file)

# add predictions
fusion_data <- add_prediction(fusion_data, model, prediction_threshold)

# Write fusions to output file -------------------------------------------------
#consitently use dec="."
write.table(fusion_data, output_file, sep=";", dec=".", row.names = FALSE)

