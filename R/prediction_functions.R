#!/usr/bin/env Rscript

# load required packages
require(tidyverse)
require(randomForest)

#' Function to replace NA  with 0 in read count columns
#' 
#' @param flat_data A data.frame with fusions. 
#' 
#' The function looks for read count columns (column names ending with "_junc" or "_span")
#' and replaces all NA's in these columns with 0. 
#' 
#' @return a data.frame with the same dimension as the input
#' 
replace_NA_in_reads <- function(flat_data){
  flat_data %>% 
    # replace NA in read count columns with zero
    mutate_at(
      .funs = ~ifelse(!is.na(.x), .x, 0),
      .vars = vars(ends_with("_junc"), ends_with("_span"))
    )
}

#' Function to change the type of some input columns to factor
#'
#' @param flat_data A data.frame with fusions.
#' @param column_names A character vector with column names that should be
#'   converted to factor.
#'
#' @return a data.frame with the same dimension as the input.
#'   
prepare_for_modeling <- function(flat_data, column_names = c("type", 
                                                             "exon_boundary", 
                                                             "frame")){
  
  flat_data %>%
    
    # convert character vecztors into factors
    mutate_at(
      .funs = as.factor,
      .vars = column_names
    ) 
}

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
  
  # get prediction probabilites
  pred_matrix <- predict(model, newdata = flat_data, type = "prob")
  pred_prob <- pred_matrix[, "positive"]
  
  # add columns
  flat_data <- flat_data %>% 
    mutate(
      prediction_prob = pred_prob,
      prediction_class = ifelse(pred_prob >= threshold, "positive", "negative")
    )
  
  return(flat_data)
}

