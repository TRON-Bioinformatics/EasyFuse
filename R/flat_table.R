
require(tidyverse)    # for analysing data in tabels
require(fs)           # to handle file paths



#' parse read counts from FASTQC file
#' 
#' @fastqc_file Path to a fastqc summary.txt file with read counts
#' 
#' 
#' 
parse_read_counts <- function(fastqc_file){
  
  readr::read_file(fastqc_file) %>% 
    str_match(pattern = "Total Sequences\\s([\\d]*)") %>% 
    magrittr::extract(, 2) %>% 
    parse_integer()
  
}

#' Reads a re-quantification file with read counts for each candidate fusion
#' 
#' @param file_path path to input file
#' 
#' Reads a requantifcation file from EasyFuse pipeline and tronsforms it into
#' a flat data table. Each line represents a unique context-sequence id. Colmns
#' with read counts have the following naming format: 
#' `<allele>_<type>` where `<allele>` is either the fusion transcript (ft), 
#' the wild-type transcript 1 (wt1) or the wild-type transcript 2 (wt2) and 
#' `<type>` is one of the follwing read cunt types:  a, b, junc, or span.
#' 
read_quantification <- function(file_path){
  
  read_delim(file_path, delim = ";", col_types = "ccciiii") %>% 
    rename(
      context_sequence_id = MD5,
      allele = Type,
      a = Reads_Gene_A, 
      b = Reads_Gene_B, 
      junc = Junc_Reads, 
      span = Span_Reads
      ) %>% 
    gather("type", "reads", a, b, junc, span) %>% 
    unite("read_type", allele, type) %>% 
    spread(read_type, reads)

  }

#' parse read counts from FASTQC file
#' 
#' @param fastqc_file Path to fastqc output file 
#' 
#' @return the total read counts as integer
#' 
parse_read_counts <- function(fastqc_file){
  content <- readr::read_file(fastqc_file) %>% 
    str_match(pattern = "Total Sequences\\s([\\d]*)") %>% 
    magrittr::extract(, 2) %>% 
    parse_integer()
}

#' A function to parse all nececary files from a EsayFuse output sample
#' 
#' @param sample_folder path to sample output folder
#' @param tools Vecotr of tool names. They should match the tools used in 
#' EasyFuse pipeline. 
#' 
#' This funciton
#' - combines potentially multiple transcripts to one entry per context sequence
#' 
#' 
flat_table <- function(sample_folder, tools = c("infusion", "starfusion", 
                                                "soapfuse", "mapsplice", 
                                                "fusioncatcher")){
  
  #-----------------------------------------------------------------------------
  # get input files from sample folder
  #-----------------------------------------------------------------------------
  detected_fusions_file <- fs::path(sample_folder, 
      "scratch/fetchdata_1tool/Detected_Fusions.csv")
  requantification_file <- fs::path(sample_folder, 
      "scratch/fetchdata_1tool/Classification.csv")
  context_seq_file <- fs::path(sample_folder, 
      "scratch/fetchdata_1tool/Context_seqs.csv")
  
  fastqc_files <- fs::dir_ls(
      fs::path(sample_folder, "scratch", "fastqc"),
      type = "directory"
    ) %>% 
    fs::path("fastqc_data.txt")

  # stop if not all files exists
  stopifnot(
    fs::file_exists(detected_fusions_file),
    fs::file_exists(requantification_file),
    fs::file_exists(context_seq_file),
    fs::file_exists(fastqc_files)
    )
  
  # assume two fastqc files 
  stopifnot(
    length(fastqc_files) == 2
  )
  

  #-----------------------------------------------------------------------------
  # read detected fusions
  #-----------------------------------------------------------------------------
  fusions_raw <- read_delim(detected_fusions_file,
                            delim = ";",
                            skip = 1,
                            col_names = c("FGID",
                                          "Fusion_Gene",
                                          "Breakpoint1",
                                          "Breakpoint2",
                                          "Junction_Reads",
                                          "Spanning_Reads",
                                          "Sample",
                                          "Tool"),
                            col_types = cols_only(
                              FGID = col_character(),
                              Fusion_Gene = col_character(),
                              Breakpoint1 = col_character(),
                              Breakpoint2 = col_character(),
                              Junction_Reads = col_integer(),
                              Spanning_Reads = col_integer(),
                              Sample = col_character(),
                              Tool = col_character()
                            ))
  
  FGID_to_annot <- fusions_raw %>%
    select(Sample, FGID, Fusion_Gene, Breakpoint1, Breakpoint2) %>%
    # TEMPORARY FIX: remove Fusion_Gene column, as Mapsplice produces duplicated
    # or multiple Gene names, inconsistend with other tools.
    select(-Fusion_Gene) %>% 
    distinct()
  
  # combine results of individual tools into a single row
  FGID_to_tool <- fusions_raw %>% 
    
    # select and rename subset of columns 
    select(Sample, FGID, tool = Tool, junc = Junction_Reads, span = Spanning_Reads) %>% 

    # add implicit missing observations of tools as NAs
    mutate(detected = TRUE) %>%
    complete(nesting(Sample, FGID), 
             tool = tools,
             fill = list("detected" = FALSE)) %>% 
    
    # combine span and junc into type colum
    gather("value_type", "value", detected, junc, span) %>% 
    unite("tool_value_type", tool, value_type) %>% 
    
    # spread in separate columns
    spread(tool_value_type, value)
   

  FGID_to_tool_count <- fusions_raw %>% 
    # add tool_count
    count(FGID) %>% 
    rename(tool_count = n) %>% 
    mutate(tool_frac = tool_count / length(tools))
  
  #-----------------------------------------------------------------------------
  # read context seq table
  #-----------------------------------------------------------------------------
  context_seq_raw <- read_delim(context_seq_file, 
                                delim = ";",
                                col_types = cols(
                                  .default = col_character(),
                                  context_sequence_bp = col_integer(),
                                  wt1_context_sequence_bp = col_integer(),
                                  wt2_context_sequence_bp = col_integer(),
                                  neo_peptide_sequence_bp = col_number()
                                )
  )
  
  # combine potentially multiple transcripts to one entry per context sequence
  context_seq <- context_seq_raw %>% 
    
    # filter out missing context_sequence
    filter(!is.na(context_sequence)) %>% 
    
    group_by(FGID, Fusion_Gene, Breakpoint1, Breakpoint2, 
             context_sequence_id, context_sequence, 
             type, exon_boundary1, exon_boundary2, exon_boundary, context_sequence_bp) %>% 
    summarize(
      n_transcripts_variants = n(),
      FTID = str_c(FTID, collapse = ";"),
      frame = str_c(frame, collapse = ";"),
      neo_peptide_sequence = str_c(neo_peptide_sequence, collapse = ";"),
      neo_peptide_sequence_bp = str_c(neo_peptide_sequence_bp, collapse = ";")
    ) %>% 
    ungroup()
    
  #-----------------------------------------------------------------------------
  # read requantification data by read calsses
  #-----------------------------------------------------------------------------
  quant_data <- read_quantification(requantification_file)
  
  #-----------------------------------------------------------------------------
  # combine all tables into a single flat table
  #-----------------------------------------------------------------------------
  flat_fusions_raw <- context_seq %>% 
    left_join(FGID_to_tool, by = "FGID") %>%
    left_join(FGID_to_tool_count, by = "FGID") %>%
    left_join(quant_data, by = c("FGID", "context_sequence_id"))
  
  # Fix mutliplie FGID and Fusion_Gene columns for same context_sequence_id
  flat_fusions <- flat_fusions_raw %>% 
    
    
    # # combine FGID and Fusion_Gene columns by ";"
    # group_by(context_sequence_id, Breakpoint1, Breakpoint2) %>% 
    # summarize(
    #   FGID = str_c(FGID, collapse = ";"),
    #   Fusion_Gene = str_c(Fusion_Gene, collapse = ";")
    # ) %>% 
    # left_join(
    #   flat_fusions_raw %>% select(-FGID, -Fusion_Gene) %>% distinct(),
    #   by = c("context_sequence_id", "Breakpoint1", "Breakpoint2")
    # )
    distinct(context_sequence_id, Breakpoint1, Breakpoint2, .keep_all = TRUE) %>% 
  
    # rearrange columns for better overview
    select(Sample, Fusion_Gene, FGID, context_sequence_id, everything())

  #-----------------------------------------------------------------------------
  # Normalize read counts to CPM
  #-----------------------------------------------------------------------------
  total_read_pairs <- map_int(fastqc_files, parse_read_counts) %>% 
    sum()
  
  #cpm = reads / total_read_pairs * 10^6  
  flat_fusions <- flat_fusions %>% 
    mutate_at(
      .vars = vars(matches("_junc$|_span$|_a$|_b$")),
      .funs = function(x){ x / total_read_pairs * 10^6}
      )

  return(flat_fusions)
  
}
