#load libraries
library(dplyr)
library(tidyr)
library(optparse)

options(stringsAsFactors = FALSE)

argument_list <- list(
        make_option(c("-i", "--input_gtf"), default="GRCh38.86_cdna_all.gtf", help="Input GTF file"),
        make_option(c("-o", "--output_csv"), default="GRCh38.86_cdna_all.csv", help="Output CSV file of transcript coordinates of the respective genome")
)

opt <- parse_args(OptionParser(option_list=argument_list))

#custome function to transform gtf
GTFtoCSV <- function(path){
	
	#read in gtf and combine info
	ensemble_gtf <- read.csv(path, header = FALSE, sep="\t", skip=5, 
		col.names = c("seqname","source","feature","start","end","score","strand","frame","attribute"))
	
	#calculate transcript information
	ensemble_transcript = ensemble_gtf %>%
		filter(feature == "transcript") %>%
		extract(attribute, into=c("transcript_id","symbol", "tsl"), 
			"transcript_id (.*); transcript_version.*gene_name (.*); gene_source.*transcript_support_level (.*)[(/;]") %>%
		mutate(gene_start = start, gene_end = end, tsl = gsub(" .*","",tsl)) %>%
		select(chr = seqname, gene_start, gene_end, strand, transcript_id, symbol, tsl) %>%
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
	ensemble_start_codon = ensemble_gtf %>%
		filter(feature == "start_codon") %>%
		extract(attribute, into="transcript_id", "transcript_id (.*); transcript_version") %>%
		select(transcript_id, start_codon = start)
	
	#calculate cds information
	ensemble_stop_codon = ensemble_gtf %>%
		filter(feature == "stop_codon") %>%
		extract(attribute, into="transcript_id", "transcript_id (.*); transcript_version") %>%
		select(transcript_id, stop_codon = end)
	
	#calculate cds information
	ensemble_cds = ensemble_gtf %>%
		filter(feature == "CDS") %>%
		extract(attribute, into="transcript_id", "transcript_id (.*); transcript_version") %>%
		group_by(transcript_id) %>%
		summarize(cds_start = min(start), cds_end = max(end)) %>%
		left_join(ensemble_start_codon, by="transcript_id") %>%
		left_join(ensemble_stop_codon, by="transcript_id") %>%
		mutate(orf = ifelse(is.na(start_codon), ifelse(is.na(stop_codon), "incomplete", "incomplete_5"), 
			ifelse(is.na(stop_codon), "incomplete_3", "complete"))) %>%
		select(transcript_id, cds_start, cds_end, orf)
		
	#combine all infos and save	
	ensemble_csv = ensemble_transcript %>%
		left_join(ensemble_cds, by="transcript_id") %>%
		left_join(rbind(ensemble_exon, ensemble_intron_plus, ensemble_intron_negative), by=c("transcript_id", "strand")) %>%
		arrange(chr, gene_start, gene_end, number, type) %>%
		group_by(chr, gene_start, gene_end, strand, transcript_id, symbol, cds_start, cds_end, tsl, orf) %>%
		summarize_all(funs(paste0(., collapse="|")))
		
	#return
	return(ensemble_csv)
}


#custome function to transform gff
GFFtoCSV <- function(path){

	#read gff file
	ensemble_gff = read.csv(path, header = FALSE, sep="\t", 
		col.names = c("seqname","source","type","start","end","score","strand","frame","attribute")) %>%
		filter(!grepl("#",seqname))

	#calculate exon information
	ensemble_exon = ensemble_gff %>%
		filter(type == "exon") %>%
		extract(attribute, into=c("transcript_id","end_phase", "phase", "number"), 
			".*transcript:(.*);Name=.*;constitutive=.;ensembl_end_phase=(.*);ensembl_phase=(.*);exon_id.*rank=(.*);version.*",convert=TRUE) %>%
		mutate(type = "E") %>%
		select(type , start, end, strand, transcript_id, end_phase, phase, number)

		
	#calculate cds information
	ensemble_cds = ensemble_gff %>%
		filter(type == "CDS") %>%
		extract(attribute, into=c("transcript_id"), 
			".*transcript:(.*);protein_id=.*",convert=TRUE) %>%
		select(cds_start = start, cds_end = end, frame, transcript_id) %>%
		left_join(ensemble_exon, by=c("transcript_id")) %>%
		filter(cds_start >= start, cds_end <= end) %>%
		select(cds_start, cds_end, frame, transcript_id, number)
	
	
	#calculate intron information
	ensemble_intron = ensemble_exon %>%
		mutate(number = number + 1) %>%
		inner_join(ensemble_exon, by=c("transcript_id", "number", "strand")) %>%
		mutate(start = ifelse(strand == "+", end.x + 1, end.y + 1),
			end = ifelse(strand == "+", start.y - 1, start.x - 1), 
			end_phase = -1, phase =-1, number = number - 0.5, type = "I") %>%
		select(type, start, end, strand, transcript_id, end_phase, phase, number)
		

	#calculate transcript information
	ensemble_csv = ensemble_gff %>%
		filter(grepl("ID=transcript", attribute)) %>%
		extract(attribute, into=c("transcript_id","symbol"), "ID=transcript:(.*);Parent=.*;Name=(.*);biotype=.*") %>%
		select(chr = seqname, gene_start = start, gene_end = end, strand, transcript_id, symbol) %>%
		left_join(rbind(ensemble_exon, ensemble_intron), by=c("transcript_id", "strand")) %>%
		left_join(ensemble_cds, by=c("transcript_id", "number")) %>%	
		arrange(chr, gene_start, gene_end, transcript_id, number) %>%
		group_by(chr, gene_start, gene_end, strand, transcript_id, symbol) %>%
		summarize_each(funs(paste0(., collapse="|")))
		
	#return
	return(ensemble_csv)	
}



#generate tables with sequences from GTF
#ensemble_csv = GTFtoCSV(path = "Mus_musculus.GRCm38.95.gtf")
#write.csv2(ensemble_csv, "GRCm38.95_cdna_all.csv", row.names=FALSE)

ensemble_csv = GTFtoCSV(path = opt$input_gtf)
write.csv2(ensemble_csv, opt$output_csv, row.names=FALSE)



#generate tables with sequences from GFF
#ensemble_csv = GFFtoCSV(path = opt$input_gff)
#write.csv2(ensemble_csv, "GRCh38.86_cdna_all_gff.csv", row.names=FALSE)

#ensemble_gtf
# SHKBP1-ENST00000600718 with incomplete gene model





