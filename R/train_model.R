# Train random forest model

# Load required packages -------------------------------------------------------

require(tidyverse)    # for analysing data in tabels
require(fs)           # to handle file paths
require(readxl)       # to read excel files tidy
require(randomForest) # for Random Forest models

# Define parameters and input files --------------------------------------------

#################### CHANGE HARD CODED PATHS HERE ##############################

PROJECT_DIR <- "/bunny/projects/CM29_RNA_Seq/fusion_gene_detection/fusion_model/"
DATA_DIR <- "/bunny/projects/CM29_RNA_Seq/fusion_gene_detection/scratch/WP4/WP4_filtered_BNT"
DATA_DIR_TRON <- "/bunny/projects/CM29_RNA_Seq/fusion_gene_detection/scratch/WP4/WP4_filtered"

################################################################################

# validation data 
validation_data_file <- path(PROJECT_DIR, "data", "2019-01-22_WP4_validation_data.csv")

# define path and prefix for output files
OUT_PREFIX = fs::path(PROJECT_DIR, "results", "Fusion_modeling_IVAC_BNT_vTEST")

# create output folders
dir.create(dirname(OUT_PREFIX), showWarnings = FALSE)

TOOLS = c("fusioncatcher", "infusion", "mapsplice", "starfusion")

sample_mapping <- tribble(
  ~BNT_id, ~sample_rep,
  "S1_fusRank_1", "0238_002_FP_MET_1TR1_RC1_S1_L001",
  "S1_fusRank_2", "0238_002_FP_MET_1TR1_RC2_S2_L001",
  "S2_fusRank_1", "0207_003_BP_MET_1TRP_SUB13_RC1_S13_L001",
  "S2_fusRank_2", "0207_003_BP_MET_1TRP_SUB13_RC2_S14_L001",
  "S3_fusRank_1", "0238_001_FP_TMA_1TR1_RC1_S9_L001",
  "S3_fusRank_2", "0238_001_FP_TMA_1TR1_RC2_S10_L001"
)

# Functions --------------------------------------------------------------------


#' Read fusion output of BNT EasyFuse pipeline as tibble
read_fusions <- function(in_file){
  read_delim(in_file, delim = ";", col_types = cols(
    .default = col_double(),
    FGID = col_character(),
    context_sequence_id = col_character(),
    FTID = col_character(),
    Fusion_Gene = col_character(),
    Breakpoint1 = col_character(),
    Breakpoint2 = col_character(),
    context_sequence_100_id = col_character(),
    type = col_character(),
    exon_starts = col_character(),
    exon_ends = col_character(),
    exon_boundary1 = col_character(),
    exon_boundary2 = col_character(),
    exon_boundary = col_character(),
    bp1_frame = col_character(),
    bp2_frame = col_character(),
    frame = col_character(),
    context_sequence = col_character(),
    neo_peptide_sequence = col_character(),
    neo_peptide_sequence_bp = col_number()
  ))
}


convert_to_factor <- function(flat_data){
  
  flat_data %>% 
    mutate(
      type = factor(type, levels = c("cis_near", "cis_inv", "cis_far", 
                                     "cis_trans", "trans", "trans_inv")),
      exon_boundary = factor(exon_boundary, levels = c("none", "3prime", 
                                                       "5prime", "both")),
      frame = factor(frame, levels = c("no_frame", "in_frame", "out_frame", 
                                       "neo_frame"))
    )
}

filter_validation <- function(flat_data){
  
  flat_data %>%
    
    # remove samples without measured validation
    filter(!is.na(validation)) %>% 
    
    # convert response variable "validation" into factor
    mutate(
      validation = factor(validation == "positive", 
                          c(FALSE, TRUE), 
                          c("negative", "positive"))
    )
}


unique_by_FGID <- function(df){
  
  # for each breakpoint take the maximum value of all context sequences
  bp_count <- df %>% 
    group_by(sample, sample_rep, FGID) %>% 
    summarize_at(.vars = vars(contains("_span"), 
                              contains("_junc"), 
                              contains("_a_"), 
                              contains("_b_"), 
                              contains("_bp_"), 
                              contains("_anch_"), 
                              ends_with("_detected"),
                              # matches("in_replicates"),
                              matches("tool_frac")), 
                 .funs = max)
  
  
  
  # consider breakpoint positive if one context sequence is positive
  bp_validation <- df %>% 
    group_by(sample, sample_rep, FGID) %>% 
    summarize(
      validation = factor("positive" %in% validation,
                          c(FALSE, TRUE),
                          c("negative", "positive")
      ) 
    )
  
  
  # combine unique columns
  bp_annot = df %>% 
    distinct(sample, sample_rep, FGID, type, exon_boundary)
  
  # combine everything to one data set
  bp_df <- bp_annot %>% 
    left_join(bp_count, by = c("sample", "sample_rep", "FGID")) %>% 
    left_join(bp_validation, by = c("sample", "sample_rep", "FGID"))
  
  # Add detection in replicates as predictive column
  bp_df <- bp_df %>% 
    add_count(sample, FGID) %>% 
    rename(in_replicates = n)
  return(bp_df)
}


# Read data --------------------------------------------------------------------

sample_tab <- sample_mapping %>% 
  mutate(sample =  str_split_fixed(sample_rep, "_RC", 2)[,1]) %>% 
  arrange(sample)

# get input files from sample folders
sample_tab <- sample_tab %>% 
  mutate(
    sample_rep_file = path(DATA_DIR, str_c(BNT_id, ".csv")),
    sample_rep_file_exists = file_exists(sample_rep_file)
  )

# read detected fusions data from sample folders
sample_tab <- sample_tab %>% 
  mutate(
    fusion_df = map(sample_rep_file, read_fusions)
  )

# combine into one large table
fusion_flat <- sample_tab %>% 
  select(sample_rep, sample, fusion_df) %>% 
  unnest(fusion_df)

# read other fusions to map context sequence IDs from 100 to 400 ---------------

sample_file <- fs::path(DATA_DIR_TRON, "samples.csv")

sample_tab_tron <- read_tsv(sample_file, col_names = c("sample_rep", "status"),
                            col_types = "cc") %>% 
  mutate(sample =  str_split_fixed(sample_rep, "_RC", 2)[,1]) %>% 
  arrange(sample) %>% 
  mutate(
    sample_rep_folder = fs::path(DATA_DIR_TRON, str_c("Sample_", sample_rep)),
    sample_rep_folder_exists = file_exists(sample_rep_folder),
    context_seq_file = fs::path(sample_rep_folder, "scratch/fetchdata_1tool/Context_seqs.csv"),
    context_seq_df = map(context_seq_file, read_delim, 
                         delim = ";",
                         col_types = cols(
                           .default = col_character(),
                           context_sequence_bp = col_integer(),
                           wt1_context_sequence_bp = col_integer(),
                           wt2_context_sequence_bp = col_integer(),
                           neo_peptide_sequence_bp = col_number()
                         )
    )
  )

context_sequence_mapping <- sample_tab_tron %>% 
  unnest(context_seq_df) %>% 
  distinct(sample, context_sequence_id, context_sequence_100_id)

# read validation data ---------------------------------------------------------

validation_data <- read_csv2(validation_data_file,
                             col_types = cols(
                               .default = col_character(),
                               X1 = col_double(),
                               product_length = col_double(),
                               cDNA = col_double(),
                               FFPE_cDNA = col_double(),
                               RT = col_double(),
                               water = col_double(),
                               Qiaxcel_Size_diff = col_double()
                             )) %>% 
  select(sample = Sample, context_sequence_id, validation = Simple_Status_auto,
         everything()) %>% 
  mutate(
    sample = sample %>% str_replace("0207_003_BP_MET_1TRP", "0207_003_BP_MET_1TRP_SUB13")
  )


# filter out NAs
validation_data_sub <- validation_data %>% 
  
  # remove potentially replicated experiments and use only selected columns
  distinct(sample, context_sequence_id, validation) %>% 
  
  # filter out NAs
  filter(!is.na(validation)) 


# add validation data to fusions -----------------------------------------------
fusion_flat <- fusion_flat %>% 
  
  # remove context sequence id, as it is only the 100 based
  select(-context_sequence_id) %>% 
  left_join(context_sequence_mapping, by = c("sample", "context_sequence_100_id")) %>% 
  
  # add validation data
  left_join(
    validation_data_sub, by = c("sample", "context_sequence_id")
  )

# add replicate count ----------------------------------------------------------
sample_FGID_rep <- fusion_flat %>% 
  distinct(sample, sample_rep, FGID) %>% 
  add_count(sample, FGID) %>% 
  rename(in_replicates = n)

fusion_flat <- fusion_flat %>% 
  left_join(sample_FGID_rep,  by =  c("sample_rep", "sample", "FGID"))


# write to output file
write_tsv(fusion_flat, str_c(OUT_PREFIX, ".fusion_flat.tsv"))


# prepare data for modelling ---------------------------------------------------

df <- fusion_flat %>% 
  convert_to_factor() %>% 
  filter_validation()

# save data
write_rds(df, str_c(OUT_PREFIX, ".df.rds"))
write_tsv(df, str_c(OUT_PREFIX, ".df.tsv"))

# consider only unique breakpoint ----------------------------------------------
bp_df <- unique_by_FGID(df)

write_rds(bp_df, str_c(OUT_PREFIX, ".bp_df.rds"))
write_tsv(bp_df, str_c(OUT_PREFIX, ".bp_df.tsv"))

# prepare feature subsets ------------------------------------------------------

annotation_features <- c("type", "exon_boundary", "in_replicates")

tool_features <- c(
  str_c(TOOLS, "_detected"),
  str_c(TOOLS, "_junc"),
  str_c(TOOLS, "_span"),
  "tool_frac"
)

alleles <- c("ft", "wt1", "wt2")

requant_junc_span <- c(
  str_c(alleles, "junc", sep = "_"),
  str_c(alleles, "span", sep = "_")
)

requant_a_b <- c(
  str_c(alleles, "a", sep = "_"),
  str_c(alleles, "b", sep = "_")
)

requant_features_flt <- str_c(c(requant_junc_span, requant_a_b), "_fltr")
requant_features_org <- str_c(c(requant_junc_span, requant_a_b), "_org")

# get predictor variables as annotated in metadata
all_predictors <- c(
  annotation_features,
  tool_features,
  requant_features_flt,
  requant_features_org
)

# get predictor variables that do not depend directely on one of the used tools

requant_breakpoint_flt <- c("exon_boundary", str_c(requant_junc_span, "_fltr"))
requant_breakpoint_org <- c("exon_boundary", str_c(requant_junc_span, "_org"))


requant_and_boundary_flt <- c("exon_boundary", requant_features_flt)
requant_and_boundary_org <- c("exon_boundary", requant_features_org)

requant_boundary_rep_flt <- c("exon_boundary", "in_replicates", requant_features_flt)
requant_boundary_rep_org <- c("exon_boundary", "in_replicates", requant_features_org)

IVAC_tools_flt <- c(requant_boundary_rep_flt, 
                    str_subset(tool_features, "starfusion|fusioncatcher"))
IVAC_tools_org <- c(requant_boundary_rep_org, 
                    str_subset(tool_features, "starfusion|fusioncatcher"))


feature_df <- tibble(
  name = c("all", 
           "requant_breakpoint_flt",
           "requant_breakpoint_org",
           "requant_and_boundary_flt",
           "requant_and_boundary_org",
           "requant_boundary_rep_flt",
           "requant_boundary_rep_org",
           "IVAC_tools_flt",
           "IVAC_tools_org"
  ),
  feature_set = list(all_predictors,
                     requant_breakpoint_flt,
                     requant_breakpoint_org,
                     requant_and_boundary_flt,
                     requant_and_boundary_org,
                     requant_boundary_rep_flt,
                     requant_boundary_rep_org,
                     IVAC_tools_flt,
                     IVAC_tools_org
  )
)

# train model ------------------------------------------------------------------


# train random forest clasifier
model_df <- feature_df %>%
  mutate(
    model = map(feature_set, ~randomForest(
      x = select(bp_df, .x), 
      y = pull(bp_df, validation), importance = TRUE),
      maxnodes = 8
    )
  )

# deploy models ----------------------------------------------------------------

# deploy all models
model_df %>% 
  mutate(
    model = walk2(model, name, ~ write_rds(.x, str_c(OUT_PREFIX, ".model.", .y, ".rds")))
  )


