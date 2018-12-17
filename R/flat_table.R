
require(tidyverse)    # for analysing data in tabels
require(fs)           # to handle file paths


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

#' A function to combine several output files into a singel table
#'
#' @param detected_fusion_file path to detected fusions file.
#' @param context_seq_file path to annotated context sequence list of detected
#'   fusions.
#' @param requantification_file path to file with requantification results.
#' @param qc_table_file path to FastQC file.
#' @param tool_state_file path to tool_state file.
#'
#'   This funciton read all detected fusions and combines it with annotations
#'   and requantificatoin data. Potentially multiple transcripts per fusion
#'   breakpoint are combined to one entry per context sequence.
#'
#' @return A data.frame with one fusion (unique context sequence per breakpoint)
#'   per row and several annotations in columns.
#' 
read_fusion <- function(detected_fusion_file, context_seq_file, requantification_file, qc_table_file, 
                        tool_state_file){
  
  #-----------------------------------------------------------------------------
  # get input files from sample folder
  #-----------------------------------------------------------------------------

  # stop if not all files exists
  stopifnot(
    fs::file_exists(detected_fusion_file),
    fs::file_exists(context_seq_file),
    fs::file_exists(requantification_file),
    fs::file_exists(qc_table_file),
    fs::file_exists(tool_state_file)
  )
  
  # parse used tools from tools_state_file
  tools <- read_lines(tool_state_file, n_max = 1) %>% 
    str_split(",") %>% 
    unlist() %>% 
    setdiff("Sample ID")
  
  message("INFO: Found the following used tools: ", str_c(tools, collapse= ",  "))

  #-----------------------------------------------------------------------------
  # read detected fusions
  #-----------------------------------------------------------------------------
  fusions_raw <- read_delim(detected_fusion_file,
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
             context_sequence_id, context_sequence_100_id, context_sequence, 
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
    
    # take only unique entries for context_sequnce and breakpoints
    distinct(context_sequence_id, Breakpoint1, Breakpoint2, .keep_all = TRUE) %>% 
  
    # rearrange columns for better overview
    select(Sample, Fusion_Gene, FGID, context_sequence_id, everything())

  #-----------------------------------------------------------------------------
  # Normalize read counts to CPM
  #-----------------------------------------------------------------------------
  total_read_pairs <- read_csv(qc_table_file, col_types = cols(
      filename = col_character(),
      read_length = col_integer(),
      suggested_trim_length = col_integer(),
      remaining_read_length = col_integer(),
      actual_trim_length = col_integer(),
      badass_bases = col_integer(),
      total_sequences = col_integer()
    )) %>% 
    pull(total_sequences) %>% 
    sum()
  
  # normlaize to counts per millions: cpm = reads / total_read_pairs * 10^6  
  flat_fusions <- flat_fusions %>% 
    mutate_at(
      .vars = vars(matches("_junc$|_span$|_a$|_b$")),
      .funs = function(x){ x / total_read_pairs * 10^6}
      )
  
  message("INFO: Combined data for ", nrow(flat_fusions), " fusions.")
  
  return(flat_fusions)
  
}
