#!/usr/bin/env Rscript

#load libraries
library("dplyr")
library("tidyr")
library("optparse")

options(stringsAsFactors = FALSE)

argument_list <- list(
  make_option(c("-i", "--input_ensembl_gtf"), default="", help="Input list of detected fusions"), 
  make_option(c("-o", "--output_csv"), default="", help="Output file of context sequences")
)
opt <- parse_args(OptionParser(option_list=argument_list))

# check mandatory arguments
if(is.na(opt$input_ensembl_gtf) | opt$input_ensembl_gtf == "") {
  print("Mandatory parameter \"input ensembl gtf\" missing. Aborting...")
}
if(is.na(opt$output_csv) | opt$output_csv == "") {
  opt$output_csv <- paste0(tools::file_path_sans_ext(opt$input_ensembl_gtf),".csv")
}


#custome function to transform gtf
GTFtoCSV <- function(path){
	
	#read in gtf and combine info
	ensemble_gtf <- read.csv(path, header = FALSE, sep="\t", skip=5,
		col.names = c("seqname","source","feature","start","end","score","strand","frame","attribute"))
	
	#calcualte transcript information
	ensemble_transcript = ensemble_gtf %>%
		filter(feature == "transcript") %>%
		extract(attribute, into=c("transcript_id","symbol"), "transcript_id (.*); transcript_version.*gene_name (.*); gene_source") %>%
		mutate(gene_start = start, gene_end = end) %>%
		select(chr = seqname, gene_start, gene_end, strand, transcript_id, symbol) %>%
		mutate(chr = as.character(chr))

	#calcualte exon information
	ensemble_exon = ensemble_gtf %>%
		filter(feature == "exon") %>%
		extract(attribute, into=c("transcript_id","number"), "transcript_id (.*); transcript_version.*exon_number (.*); gene_name", convert=TRUE) %>%
		mutate(type = "exon") %>%
		select(start, end, number, strand, transcript_id, type)
	
	#calculate intron information
	ensemble_intron_plus = ensemble_exon %>%
		filter(strand == "+") %>%
		mutate(number = number + 1) %>%
		inner_join(ensemble_exon, by=c("transcript_id", "number", "strand")) %>%
		transmute(start = end.x + 1, end = start.y - 1, number = number - 0.5, strand = strand, transcript_id = transcript_id, type = "intron")
	
	ensemble_intron_negative = ensemble_exon %>%
		filter(strand == "-") %>%
		mutate(number = number + 1) %>%
		inner_join(ensemble_exon, by=c("transcript_id", "number", "strand")) %>%
		transmute(start = end.y + 1, end = start.x - 1, number = number - 0.5, strand = strand, transcript_id = transcript_id, type = "intron")
		
	#calculate cds information
	ensemble_cds = ensemble_gtf %>%
		filter(feature == "CDS") %>%
		extract(attribute, into="transcript_id", "transcript_id (.*); transcript_version") %>%
		group_by(transcript_id) %>%
		summarize(cds_start = min(start), cds_end = max(end))
		
	#combine all infos and save	
	ensemble_csv = ensemble_transcript %>%
		left_join(ensemble_cds, by="transcript_id") %>%
		left_join(rbind(ensemble_exon, ensemble_intron_plus, ensemble_intron_negative), by=c("transcript_id", "strand")) %>%
		arrange(chr, gene_start, gene_end, number, type) %>%
		group_by(chr, gene_start, gene_end, strand, transcript_id, symbol, cds_start, cds_end) %>%
		summarize_each(funs(paste0(., collapse="|")))
		
	#return
	return(ensemble_csv)
}

#generate tables with sequences
ensemble_csv = GTFtoCSV(path = opt$input_ensembl_gtf)
write.csv2(ensemble_csv, opt$output_csv, row.names=FALSE)
