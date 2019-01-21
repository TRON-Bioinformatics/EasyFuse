#load libraries
#print(.libPaths())
# urla: R doesn't look in the correct path for the icui libs when run from anaconda.
#       run "export LD_LIBRARY_PATH=/path/to/lib:$LD_LIBRARY_PATH" in the exe terminal to circumvent the error
#       e.g. export LD_LIBRARY_PATH=/data/urla_progs/anaconda2/lib:$LD_LIBRARY_PATH
library("Biostrings")
library("GenomicRanges")
library("dplyr")
library("tidyr")
library("digest")
library("BSgenome")
library("optparse")

options(stringsAsFactors = FALSE)

argument_list <- list(
	make_option(c("-i", "--input_detected_fusions"), default="Detected_Fusions.csv", help="Input list of detected fusions"), 
	make_option(c("-f", "--fasta_genome_dir"), default="", help="Directory containing fasta sequences (one file per chr/contig) of the respective genome"), 
	make_option(c("-e", "--ensembl_csv"), default="GRCh38.86_cdna_all.csv", help="Gtf derived csv file of transcript coordinates of the respective genome"), 
	make_option(c("-o", "--output"), default="context_seqs.csv", help="Output file of context sequences"), 
	make_option("--cis_near_distance", type="integer", default=1000000, help="Distance to call fusion cis near or far"), 
	make_option("--genomic_seq_len", type="integer", default=1000, help="no idea"), 
	make_option("--context_seq_len", type="integer", default=100, help="Length of the context sequence")
)
opt <- parse_args(OptionParser(option_list=argument_list))

# check mandatory arguments
if(is.na(opt$input_detected_fusions) | opt$input_detected_fusions == "") {
	print("Mandatory parameter \"detected fusions file\" missing. Aborting...")
}
if(is.na(opt$fasta_genome_dir) | opt$fasta_genome_dir == "") {
	print("Mandatory parameter \"fasta genome directory\" missing. Aborting...")
}
if(is.na(opt$ensembl_csv) | opt$ensembl_csv == "") {
	print("Mandatory parameter \"ensembl csv file\" missing. Aborting...")
}
if(is.na(opt$output) | opt$output == "") {
	print("Mandatory parameter \"output file\" missing. Aborting...")
}

in_file = opt$input_detected_fusions
out_file = opt$output
fa_path = opt$fasta_genome_dir #fa_path = paste0(getwd(), "/hg38_chroms/")
ensemble_csv_file = opt$ensembl_csv
cis_near_distance = opt$cis_near_distance
genomic_seq_len = opt$genomic_seq_len
context_seq_len = opt$context_seq_len


# select optimal type of matching transcript
select_type <- function(var1){
	if(any(grepl("boundary", var1))) {
		return("boundary")
	} else if(any(grepl("exon", var1))) {
		return("exon")
	} else if(any(grepl("intron", var1))) {
		return("intron")
	} else {
		return("intergenic")
	}
}



#custom function read all genomic sequences - using granges
read_genome_seq <- function(chr, start, end, strand){

	#list genome fasta files
	chr_names = gsub(".fa", "", list.files(fa_path, pattern = "*.fa$", full.names = FALSE))
	chr_files = list.files(fa_path, pattern = "*.fa$", full.names = TRUE)
	chr_dat = data.frame(chr = chr_names, chr_files)
	
	# urla: I don't like dplyr pipes in code as it makes it very hard to read and understand!
	# dawb: Personally I think it is much easier to read, furthermore execution of tidyr/dplyr 
	# is much faster than base R, I would not recommend to change it

	tmp_data = data.frame(chr, start, end, strand) %>%
    left_join(chr_dat, by="chr")
    
	#empty dataframe for results
	seq = c()
	
	#GenomicRanges
	for (i in unique(tmp_data$chr)){
		i_ranges = makeGRangesFromDataFrame(filter(tmp_data, chr == i))
			
		print(paste0("->load chr: ", i))
		chr_seq = tryCatch(readDNAStringSet(filter(chr_dat, chr == i)$chr_files), error = function(err) return(""))
	
		if (chr_seq != ""){
			names(chr_seq) = i #rename chr
			i_seq = as.character(getSeq(chr_seq, i_ranges)) #extract sequences
		}else{
			i_seq = rep("", length(i_ranges))
		}
		
		seq = c(seq, i_seq)
		
	}#for i

	return(seq)
}



#determine relative distance of exon end to tss by adding up sizes - not in use currently
calculate_rel_start <- function(FGID, transcript_id, number, size){
	tmp_data = data.frame(FGID, transcript_id, number, size)
	
	rel_end = c()
	for (i in 1:nrow(tmp_data)){
		tmp = tmp_data %>%
			filter(FGID == tmp_data$FGID[i], transcript_id == tmp_data$transcript_id[i], number < tmp_data$number[i]) #filter all exons 1:n for current transcript and current exon n
		
		rel_end = c(rel_end, sum(tmp$size)) #sum up exon sizes
		
	}

	return(rel_end)
}



#determine the frame acording to distance in cds
calculate_frame <- function(bp_type, side, cds_start_rel, cds_end_rel, bp_rel){
	tmp_data = data.frame(bp_type, side, cds_start_rel, cds_end_rel, bp_rel, frame = "NA") %>%
		mutate(bp_rel = ifelse(side == 2, bp_rel - 1, bp_rel))
    
	for (i in 1:nrow(tmp_data)){
	
		if(tmp_data$bp_type[i] == "intergenic"){
			tmp_data$frame[i] = "intergenic"
		}else if(tmp_data$bp_type[i] == "intron"){
			tmp_data$frame[i] = "intron"
		}else{
			if(is.na(tmp_data$cds_start_rel[i])){
				tmp_data$frame[i] = "non-coding"
			}else{
				if(tmp_data$cds_start_rel[i] > tmp_data$bp_rel[i]){
					tmp_data$frame[i] = "5UTR"
				}else{
					if(tmp_data$cds_end_rel[i] < tmp_data$bp_rel[i]){
						tmp_data$frame[i] = "3UTR"
					}else{
						tmp_data$frame[i] = (tmp_data$bp_rel[i] - tmp_data$cds_start_rel[i]) %% 3
					}
				}
			}
		}
	}
	
	return(tmp_data$frame)
}

	

#custom function to determine frame of fusion according to frames of WT transcripts
determine_frame <- function(S1_frame, S2_frame){
	tmp_data = data.frame(S1_frame, S2_frame) %>%
		mutate(frame = ifelse(S1_frame %in% c("0", "1", "2", "3"), ifelse(S2_frame %in% c("0", "1", "2", "3"),
			ifelse(S1_frame == S2_frame, "in_frame", "out_frame"), "neo_frame"), "no_frame"))
	  
	return(tmp_data$frame)
}



#custom function to determine type of bp_position
determine_type <- function(S1_chr, S1_pos, S1_strand, S2_chr, S2_pos, S2_strand){
	tmp_data = data.frame(S1_chr, S1_pos, S1_strand, S2_chr, S2_pos, S2_strand) %>%
		mutate(type = ifelse(S1_chr != S2_chr, 
			ifelse(S1_strand != S2_strand, "trans_inv", "trans"), ifelse(S1_strand != S2_strand, "cis_inv", 
			ifelse(S1_strand == "+", ifelse(S2_pos - S1_pos < 0, "cis_trans", ifelse(S2_pos - S1_pos > cis_near_distance, "cis_far", "cis_near")), 
			ifelse(S1_pos - S2_pos < 0, "cis_trans", ifelse(S1_pos - S2_pos > cis_near_distance, "cis_far", "cis_near"))))))
  
	return(tmp_data$type)
}



#custom function to translate into peptide
pep_translate <- function(cds_seq){
	pep = c()
	for (i in 1:length(cds_seq)){
		pep = c(pep, tryCatch({sub("\\*.*", "", as.character(translate(DNAString(cds_seq[i]), genetic.code = GENETIC_CODE, if.fuzzy.codon = "error")))},
		               error=function(err) {""}))
	}
	
	return(pep)
}



#custom function to generate context_seq_id
make_context_id <- function(context_sequence){
	id = c()
	for (i in 1:length(context_sequence)){
		id = c(id, digest(context_sequence[i], "xxhash64", serialize = FALSE)) #no vectorized function
	}
	
	return(id)
}



########################################################
###### procesing of detected fusions starts here #######
########################################################

tmp_time = Sys.time()

#read in and split breakpoints
print(paste0("->read detected fusion from: ", in_file))
breakpoints = read.csv(in_file, sep = ";", dec = ".", na.strings = c("na")) %>%
	distinct(FGID, Fusion_Gene, Breakpoint1, Breakpoint2) %>%
	gather(side, Breakpoint, -FGID, -Fusion_Gene) %>%
	mutate(side = as.numeric(gsub("Breakpoint", "", side))) %>%
	separate(Breakpoint, sep = ":", into = c("chr", "pos", "strand"), convert = TRUE) %>%
	mutate(chr = gsub("chr", "", chr), strand = ifelse(strand == ".", "+|-", strand)) %>%
	unnest(strand = strsplit(strand, "\\|")) %>%
	mutate(element = ifelse((strand == "+" & side == 1) | (strand == "-" & side == 2), 1, 2)) #element 1 keep left side, element 2 keep right side genomic



	
#read in ensemble data
print(paste0("->read in transcript information from: ", ensemble_csv_file))
ensemble_csv = breakpoints %>%
	transmute(chr = chr,
		gene_start = pos - genomic_seq_len, gene_end = pos + genomic_seq_len,
		strand = strand, transcript_id = paste0("intergenic_", chr, ":", pos, ":", strand), symbol = "intergenic",
		cds_start = NA, cds_end = NA, 
		start = ifelse(element == 1, paste0(pos - genomic_seq_len, "|", pos + 1), paste0(pos - genomic_seq_len, "|", pos)),
		end = ifelse(element == 1, paste0(pos, "|", pos + genomic_seq_len), paste0(pos - 1, "|", pos + genomic_seq_len)),
		number = ifelse(strand == "+", "1|2", "2|1"), type = "intergenic|intergenic") %>% #create 2000bp "intergenic genes" for each breakpoint with two 1000bp elements, as alternative when no gene is found
	rbind(read.csv2(ensemble_csv_file))



#identify overlapping elements on genome coordinates
print("->identify overlapping elements")
genomic_overlap = breakpoints %>% 
  #filter(FGID == "C1ORF100_chr1:244363704:+_SNAP47_chr1:227735499:+") %>%
	inner_join(ensemble_csv, by = c("chr","strand")) %>%   # <- huge memory usage ???
	filter(pos >= gene_start, pos <= gene_end) %>% #filter overlapping genes
	unnest(start2 = strsplit(start, "\\|"), end2 = strsplit(end, "\\|"), bp_element = strsplit(number, "\\|"), bp_type = strsplit(type, "\\|")) %>% #unnest bp_elements of overlapping genes
	mutate(start2 = as.numeric(start2), end2 = as.numeric(end2), bp_element = as.numeric(bp_element)) %>%
	filter(pos >= start2, pos <= end2) %>% #filter overlapping bp_elements
	mutate(bp_type = ifelse((element == 1 & bp_type == "exon" & pos == end2)|
		(element == 2 & bp_type == "exon" & pos == start2), "boundary", bp_type)) %>% #mark bp_elements with matching exon boundaries
	group_by(FGID, Fusion_Gene, side, chr, pos, strand) %>%
	filter(bp_type == select_type(bp_type)) %>% #filter best matching transcripts per breakpoint
	arrange() %>%
  ungroup() %>%
  distinct(FGID, Fusion_Gene, side, element, chr, pos, strand, bp_type, bp_element, transcript_id, symbol, 
		gene_start, gene_end, cds_start, cds_end, start, end, number, type) %>%
	unnest(start = strsplit(start, "\\|"), end = strsplit(end, "\\|"), number = strsplit(number, "\\|"), type = strsplit(type, "\\|")) %>%
	mutate(start = as.numeric(start), end = as.numeric(end), number = as.numeric(number)) %>%
	filter(type != "intron" | number == bp_element) %>% #keep only exons and introns and intergenic elements with breakpoints
	mutate(size = end - start + 1)



#load sequence for all elements
print(paste0("->generate transcript sequences from: ", fa_path))
transcript_seq = genomic_overlap %>%
  distinct(FGID, transcript_id, chr, start, end, strand, number) %>%
  arrange(chr, transcript_id, number) %>% #reorder to reduce loading and ensure correct assembly of exon sequences to transcript
	mutate(seq = read_genome_seq(chr, start, end, strand)) %>% #determine sequence of all elements
	group_by(FGID, transcript_id) %>%
	summarize(seq = paste0(seq, collapse = ""), exon_nr = n(), exon_starts = paste0(start, collapse="|"), exon_ends = paste0(end, collapse="|")) %>%
	ungroup()



#calculate relative positions of elements/exons on transcript
print("->calculate relative positions on transcript")
transcript_overlap_elements = genomic_overlap %>%	
	full_join(select(genomic_overlap, FGID, transcript_id, number2 = number, size2 = size), by = c("FGID", "transcript_id")) %>% # combine all exons per transcript
	filter(number >= number2) %>%
	group_by(FGID, transcript_id, number, size) %>%
	summarize(end_rel = sum(size2)) %>%
  ungroup() %>%
  mutate(start_rel = end_rel - size)



	
#calculate relative positions of cds
transcript_cds = genomic_overlap %>%
	filter((cds_start >= start & cds_start <= end)|(cds_end >= start & cds_end <= end)) %>% #filter elements with cds_start or cds_end
	inner_join(transcript_overlap_elements, by = c("FGID", "transcript_id", "number", "size")) %>%
	mutate(cds_start_rel = ifelse(strand == "+", cds_start - start + start_rel + 1, end - cds_end + start_rel + 1),
	   cds_end_rel = ifelse(strand == "+", cds_end - start + start_rel, end - cds_start + start_rel)) %>% #calculate 
	group_by(FGID, transcript_id) %>%
	summarize(cds_start_rel = max(cds_start_rel), cds_end_rel = min(cds_end_rel)) %>%
  ungroup()



#calculate relative positions of bp and frame and combine with transcript sequence 
transcript_bp = genomic_overlap %>%
	filter(bp_element == number) %>%
	inner_join(transcript_overlap_elements, by=c("FGID", "transcript_id", "number", "size")) %>%
	left_join(transcript_cds, by=c("FGID", "transcript_id")) %>%
	mutate(pos_rel = ifelse(strand == "+", pos - start + start_rel + 1, end - pos + start_rel + 1)) %>% #calculate relative distance to bp account for strand
	mutate(frame = calculate_frame(bp_type, side, cds_start_rel, cds_end_rel, pos_rel)) %>% #determine frame
	inner_join(transcript_seq, by = c("FGID", "transcript_id"))
	


#filter information for side 1 of the breakpoint
overlap_side1 = filter(transcript_bp, side == 1) #breakpoint1
colnames(overlap_side1)[3:ncol(overlap_side1)] = paste0("S1_", colnames(overlap_side1)[3:ncol(overlap_side1)])

#filter information for side 2 of the breakpoint
overlap_side2 = filter(transcript_bp, side == 2) #breakpoint2
colnames(overlap_side2)[3:ncol(overlap_side2)] = paste0("S2_", colnames(overlap_side2)[3:ncol(overlap_side2)])



#combine both sides and extract relevant information
print("->combine all information")
overlap = overlap_side1 %>%
	full_join(overlap_side2, by = c("FGID", "Fusion_Gene")) %>%
	distinct() %>%
	transmute(FGID = FGID, Fusion_Gene, 
		Breakpoint1 = paste0(S1_chr, ":", S1_pos, ":", S1_strand), Breakpoint2 = paste0(S2_chr, ":", S2_pos, ":", S2_strand), #determine breakpoints
		FTID = paste0(S1_symbol, "_", S1_chr, ":", S1_pos, ":", S1_strand,"_", S1_transcript_id,  "_", S2_symbol, "_", S2_chr, ":", S2_pos, ":", S2_strand, "_", S2_transcript_id), #old style FTID
		#FTID = paste0(S1_symbol, "_", S1_transcript_id, "_", S1_chr, ":", S1_pos, ":", S1_strand, "_", S2_symbol, "_", S2_transcript_id, "_", S2_chr, ":", S2_pos, ":", S2_strand), #calculate fusion transcript ID
		type = determine_type(S1_chr, S1_pos, S1_strand, S2_chr, S2_pos, S2_strand), #determine how breakpoints match	
		exon_boundary1 = S1_bp_type, exon_boundary2 = S2_bp_type, 
		exon_boundary = ifelse(S1_bp_type == "boundary", ifelse(S2_bp_type == "boundary", "both", "5prime"), ifelse(S2_bp_type == "boundary", "3prime", "none")), #determine if bp matches with known exon_boundaries
		bp1_frame = S1_frame, bp2_frame = S2_frame,	frame = determine_frame(S1_frame, S2_frame),
		#wt1_sequence = S1_seq, wt1_sequence_bp = S1_pos_rel, wt2_sequence = S2_seq, wt2_sequence_bp = S2_pos_rel, #determine wildtype transcripts and breakpoint positions
		exon_nr = S1_exon_nr + S2_exon_nr, exon_starts = paste0(S1_exon_starts, "|", S2_exon_starts), exon_ends = paste0(S1_exon_ends, "|", S2_exon_ends), #combine exon positions
		transcript_sequence = paste0(substr(S1_seq, 1, S1_pos_rel), substr(S2_seq, S2_pos_rel, nchar(S2_seq))), transcript_sequence_bp = S1_pos_rel, #combine sequence for fusion transcript
		wt1_context_sequence = substr(S1_seq, S1_pos_rel - context_seq_len + 1, S1_pos_rel + context_seq_len), wt1_context_sequence_bp = pmin(S1_pos_rel, context_seq_len),
		wt2_context_sequence = substr(S2_seq, S2_pos_rel - context_seq_len, S2_pos_rel + context_seq_len - 1), wt2_context_sequence_bp = pmin(S2_pos_rel, context_seq_len),
		context_sequence_100 = paste0(substr(S1_seq, S1_pos_rel - 100, S1_pos_rel), substr(S2_seq, S2_pos_rel, S2_pos_rel + 100 - 1)), #determine context_seq_100 sequence with fixed length
		context_sequence_100_bp = pmin(S1_pos_rel, 100 + 1),
		context_sequence = paste0(substr(S1_seq, S1_pos_rel - context_seq_len, S1_pos_rel), substr(S2_seq, S2_pos_rel, S2_pos_rel + context_seq_len - 1)), #determine context_seq sequence
		context_sequence_bp = pmin(S1_pos_rel, context_seq_len + 1),
		cds_sequence = ifelse(S1_frame %in% c(0, 1, 2, 3), paste0(substr(S1_seq, S1_cds_start_rel, S1_pos_rel), substr(S2_seq, S2_pos_rel, nchar(S2_seq))), ""), #determine coding sequence 
		cds_sequence_bp = ifelse(S1_frame %in% c(0, 1, 2, 3), S1_pos_rel - S1_cds_start_rel, NA)) %>%
	mutate(context_sequence_100_id = make_context_id(context_sequence_100), context_sequence_id = make_context_id(context_sequence), #generate context seq ids
		peptide_sequence = pep_translate(cds_sequence), peptide_sequence_bp = round(cds_sequence_bp / 3,1)) %>% #translate cds to peptide
	mutate(neo_peptide_sequence = ifelse(frame == "in_frame", 
		substr(peptide_sequence, round(peptide_sequence_bp - 13.5, 0), round(peptide_sequence_bp + 13.5, 0)),
		substr(peptide_sequence, round(peptide_sequence_bp - 13.5, 0), nchar(peptide_sequence))), #determine sequence potentially immunogenic
		neo_peptide_sequence_bp = peptide_sequence_bp - round(peptide_sequence_bp - 13.5, 0))



# select relevant columns for context_seq.csv
output = overlap %>%
	select(FGID, Fusion_Gene, Breakpoint1, Breakpoint2, FTID, context_sequence_id, context_sequence_100_id,
		type, exon_nr, exon_starts, exon_ends, exon_boundary1, exon_boundary2, exon_boundary, bp1_frame, bp2_frame, frame, 
		context_sequence, context_sequence_bp, neo_peptide_sequence, neo_peptide_sequence_bp)

write.csv2(output, out_file, row.names = FALSE, quote = FALSE)



# urla: write context seqs to a fasta file with the Biostrings library function writeXStringSet and write context seq based bed file for classification with R base functions
# urla: see comment in original "custom_transcriptome.py on the "problem" with bed file generation
# dawb: further fasta files are generated for plausibility checks

# dataframe of sequences for fasta file
print("->generate fasta data")
context_seq_data = overlap %>%
	transmute(name = paste0(FTID, "_", context_sequence_id), ft_context_sequence = paste0(context_sequence, "_", context_sequence_bp), 
	  wt1_context_sequence = paste0(wt1_context_sequence, "_", wt1_context_sequence_bp), wt2_context_sequence = paste0(wt2_context_sequence,"_", wt2_context_sequence_bp)) %>%
	gather(type, sequence, -name) %>%
    separate(sequence, into=c("sequence","bp"), sep="_") %>%
	mutate(name = paste0(name, "_", bp, "_", gsub("_context_sequence","",type))) %>%
	distinct(name, sequence) %>%
	arrange(name) 

seqSet <- DNAStringSet(context_seq_data$sequence)
names(seqSet) <- context_seq_data$name
writeXStringSet(seqSet, paste0(out_file, ".fasta"), append = F, format = "fasta", width = 60)
write.table(paste(length(seqSet), sum(width(seqSet)), sep = "\n"), paste(out_file, ".fasta.info", sep = ""), row.names = F, col.names = F, quote = F)

transcrtipt_seq_data = overlap %>%
	transmute(name = paste0(FTID,"_", context_sequence_id), ft_sequence = paste0(transcript_sequence,"_",transcript_sequence_bp),
		cds_sequence = paste0(cds_sequence,"_",cds_sequence_bp)) %>%
	gather(type, sequence, -name) %>%
	separate(sequence, into=c("sequence","bp"), sep="_") %>%
	mutate(name = paste0(name, "_", gsub("_sequence","",type),":",bp)) %>%
	distinct(name, sequence) %>%
	arrange(name)


transcriptseqSet <- DNAStringSet(transcrtipt_seq_data$sequence)
names(transcriptseqSet) <- transcrtipt_seq_data$name
writeXStringSet(transcriptseqSet, paste0(out_file, "_transcript.fasta"), append = F, format = "fasta", width = 60)


peptide_seq_data = overlap %>%
	transmute(name = paste0(FTID,"_", context_sequence_id, "_", peptide_sequence_bp, "_", frame), 
		sequence = peptide_sequence) %>%
	distinct(name, sequence)

peptideseqSet <- AAStringSet(peptide_seq_data$sequence)
names(peptideseqSet) <- peptide_seq_data$name
writeXStringSet(peptideseqSet, paste0(out_file, "_peptide.fasta"), append = F, format = "fasta", width = 60)



# dataframe with positions for bed file
print("->generate bed data")
bed_data = overlap %>%
	extract(FGID, into=c("Gene1","Strand1","Gene2","Strand2"), regex="^(.*)_.*([//+|//-])_(.*)_.*([//+|//-])", remove=FALSE) %>% #extract gene and strand
	group_by(ID = paste0(FGID,"_",context_sequence_id)) %>% #remove redundancy per ID, keep only first
	summarize(Gene = first(paste0(Gene1,"|",Gene2)), Strand = first(paste0(Strand1,"|",Strand2)), 
		ft = first(paste0("0,", context_sequence_bp,"|", context_sequence_bp,",", nchar(context_sequence))),
		wt1 = first(paste0("0,", wt1_context_sequence_bp,"|", wt1_context_sequence_bp,",", nchar(wt1_context_sequence))),
		wt2 = first(paste0("0,", wt2_context_sequence_bp,"|", wt2_context_sequence_bp,",", nchar(wt2_context_sequence)))) %>% #calculate positions
	gather(type, Position, -ID, -Gene, -Strand) %>%
	arrange(ID, type) %>%
	unnest(Strand = strsplit(Strand, "\\|"), Gene = strsplit(Gene, "\\|"), Position = strsplit(Position, "\\|")) %>%
	separate(Position, into=c("startA", "endA"), sep=",", convert=TRUE) %>%
	transmute(ID = paste0(ID, "_", type), startA = startA, endA = endA, gene = Gene, emptyOnly = "", strand = Strand,
		startB = startA, endB = endA, zeroOnly=0, oneOnly=1, startC = ifelse(startA == 0, endA, endA-startA), endC = startA) 

write.table(bed_data, paste(out_file, ".bed", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)
print(paste0("->finished, overall processing time: ", Sys.time() - tmp_time))
