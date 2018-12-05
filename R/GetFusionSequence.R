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
library("bindrcpp")
library("BSgenome")
library("optparse")

options(stringsAsFactors = FALSE)

argument_list <- list(
  make_option(c("-i", "--input_detected_fusions"), default="", help="Input list of detected fusions"), 
  make_option(c("-f", "--fasta_genome_dir"), default="", help="Directory containing fasta sequences (one file per chr/contig) of the respective genome"), 
  make_option(c("-e", "--ensembl_csv"), default="", help="Gtf derived csv file of transcript coordinates of the respective genome"), 
  make_option(c("-o", "--output"), default="", help="Output file of context sequences"), 
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
fa_path = opt$fasta_genome_dir
ensemble_csv_file = opt$ensembl_csv
cis_near_distance = opt$cis_near_distance
genomic_seq_len = opt$genomic_seq_len
context_seq_len = opt$context_seq_len

#win test
#in_file = "C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\Scripting\\testGetFusionSequence\\data\\Detected_Fusions.csv"
#out_file = "C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\Scripting\\testGetFusionSequence\\data\\context_seqs.csv"
#fa_path = "C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\Scripting\\testGetFusionSequence\\data\\chromosomes_PAG"
#ensemble_csv_file = "C:\\Users\\Urs.Lahrmann\\Documents\\FusionPrediction\\Scripting\\testGetFusionSequence\\data\\Homo_sapiens.GRCh38.86.csv"
#cis_near_distance = 1000000
#genomic_seq_len = 1000
#context_seq_len = 100

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
read_genome_seq2 <- function(chr, start, end, strand){
	tmp_data = data.frame(chr, start, end, strand)
		
	# #for testing
	# tmp_data = genomic_overlap %>%
		# group_by(FGID, transcript_id, chr, start, end, strand, number) %>%
		# summarize() %>%
		# ungroup() %>%
		# arrange(chr, transcript_id, number) %>% #reorder to reduce loading and ensure correct esembly of exon sequences to transcript
		# select(chr, start, end, strand)
		
	#list genome fasta files
	chr_names = gsub(".fa", "", list.files(fa_path, pattern = "*.fa$", full.names = FALSE))
	chr_files = list.files(fa_path, pattern = "*.fa$", full.names = TRUE)
	# urla: I don't like dplyr pipes in code as it makes it very hard to read and understand!
	# urla_c: is the following equivalent to
	chr_dat = data.frame(chr = chr_names, chr_files)
	chr_dat = chr_dat[chr_dat$chr %in% tmp_data$chr, ]
	# this here?
#	chr_dat = data.frame(chr = chr_names, chr_files) %>%
#		semi_join(tmp_data, by = "chr") #filter for relevant chr
	
	#empty dataframe for results
	seq = c()
	
	#GenomicRanges
	for (i in 1:nrow(chr_dat)){
#		i_ranges = tmp_data %>%
#			filter(chr == chr_dat$chr[i]) %>%#filter data for current chr
#			makeGRangesFromDataFrame()
		i_ranges = makeGRangesFromDataFrame(tmp_data[tmp_data$chr == chr_dat$chr[i], ])
			
		print(paste0("->load chr: ", chr_dat$chr[i]))
		chr_seq = tryCatch(readDNAStringSet(chr_dat$chr_files[i]), error = function(e) return(""))
	
		if (chr_seq != ""){
			names(chr_seq) = chr_dat$chr[i] #rename chr
			i_seq = as.character(getSeq(chr_seq, i_ranges)) #extract sequences
		}else{
			i_seq = ""
		}
		
		if (i == 1){
			seq = i_seq
		}else{
			seq = c(seq, i_seq)
		}
	}#for i
	
	return(seq)
}


#determine relative exon end to tss by adding up sizes - not in use currently
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
		pep = c(pep, sub("\\*.*", "", as.character(translate(DNAString(cds_seq[i]), genetic.code = GENETIC_CODE, if.fuzzy.codon = "error"))))
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
	group_by(FGID, Fusion_Gene, Breakpoint1, Breakpoint2) %>%
	summarize() %>%
	gather(side, Breakpoint, -FGID, -Fusion_Gene) %>%
	mutate(side = as.numeric(gsub("Breakpoint", "", side))) %>%
	separate(Breakpoint, sep = ":", into = c("chr", "pos", "strand"), convert = TRUE) %>%
	mutate(chr = gsub("chr", "", chr), strand = ifelse(strand == ".", "+|-", strand)) %>%
	unnest(strand = strsplit(strand, "\\|")) %>%
	mutate(element = ifelse((strand == "+" & side == 1) | (strand == "-" & side == 2), 1, 2)) %>% #element 1 keep left side, element 2 keep right side genomic
	ungroup()


	
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


	
#identify overlapping elements on genome 
print("->identify overlapping elements")
genomic_overlap = breakpoints %>% 
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
	group_by(FGID, Fusion_Gene, side, element, chr, pos, strand, bp_type, bp_element, transcript_id, symbol, 
		gene_start, gene_end, cds_start, cds_end, start, end, number, type) %>%
	summarize() %>%
	ungroup() %>%
	unnest(start = strsplit(start, "\\|"), end = strsplit(end, "\\|"), number = strsplit(number, "\\|"), type = strsplit(type, "\\|")) %>%
	mutate(start = as.numeric(start), end = as.numeric(end), number = as.numeric(number)) %>%
	filter(type != "intron" | number == bp_element) %>% #keep only exons and introns and intergenic elements with breakpoints
	mutate(size = end - start + 1)



#load sequence for all elements
print(paste0("->generate transcript sequences from: ", fa_path))
transcript_seq = genomic_overlap %>%
	group_by(FGID, transcript_id, chr, start, end, strand, number) %>%
	summarize() %>%
	ungroup() %>%
	arrange(chr, transcript_id, number) %>% #reorder to reduce loading and ensure correct esembly of exon sequences to transcript
	mutate(seq = read_genome_seq2(chr, start, end, strand)) %>% #determine sequence of all elements
	group_by(FGID, transcript_id) %>%
	summarize(seq = paste0(seq, collapse = "")) %>%
	ungroup()
  


#calculate relative positions of elements/exons on transcript
print("->calculate relative positions on transcript")
transcript_overlap_elements = genomic_overlap %>%	
	full_join(select(genomic_overlap, FGID, transcript_id, number2 = number, size2 = size), by = c("FGID", "transcript_id")) %>% # combine all exons per trasncript
	filter(number >= number2) %>%
	group_by(FGID, transcript_id, number, size) %>%
	summarize(end_rel = sum(size2)) %>%
	mutate(start_rel = end_rel - size) %>%
	ungroup()


	
#calculate relative positions of cds
transcript_cds = genomic_overlap %>%
                filter((cds_start >= start & cds_start <= end)|(cds_end >= start & cds_end <= end)) %>% #filter elements with cds_start or cds_end
                inner_join(transcript_overlap_elements, by = c("FGID", "transcript_id", "number", "size")) %>%
                mutate(cds_start_rel = ifelse(strand == "+", cds_start - start + start_rel + 1, end - cds_end + start_rel + 1),
                               #cds_end_rel = ifelse(strand == "-", cds_start - start + start_rel, end - cds_end + start_rel)) %>% #error corrected 
                               cds_end_rel = ifelse(strand == "+", cds_end - start + start_rel, end - cds_start + start_rel)) %>% #calculate 
                group_by(FGID, transcript_id) %>%
                summarize(cds_start_rel = max(cds_start_rel), cds_end_rel = min(cds_end_rel))




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

#filter information for side 1 of the breakpoint
overlap_side2 = filter(transcript_bp, side == 2) #breakpoint2
colnames(overlap_side2)[3:ncol(overlap_side2)] = paste0("S2_", colnames(overlap_side2)[3:ncol(overlap_side2)])


#combine both sides and extract relevant information
print("->combine all information")
overlap = overlap_side1 %>%
	full_join(overlap_side2, by = c("FGID", "Fusion_Gene")) %>%
	distinct() %>%
	mutate(FGID = FGID, Fusion_Gene, 
		Breakpoint1 = paste0(S1_chr, ":", S1_pos, ":", S1_strand), Breakpoint2 = paste0(S2_chr, ":", S2_pos, ":", S2_strand), #determine breakpoints
		FTID = paste0(S1_symbol, "_", S1_chr, ":", S1_pos, ":", S1_strand,"_", S1_transcript_id,  "_", S2_symbol, "_", S2_chr, ":", S2_pos, ":", S2_strand, "_", S2_transcript_id), #old style FTID
		#FTID = paste0(S1_symbol, "_", S1_transcript_id, "_", S1_chr, ":", S1_pos, ":", S1_strand, "_", S2_symbol, "_", S2_transcript_id, "_", S2_chr, ":", S2_pos, ":", S2_strand), #calculate fusion transcript ID
		type = determine_type(S1_chr, S1_pos, S1_strand, S2_chr, S2_pos, S2_strand), #determine how breakpoints match	
		exon_boundary1 = S1_bp_type, exon_boundary2 = S2_bp_type, 
		exon_boundary = ifelse(S1_bp_type == "boundary", ifelse(S2_bp_type == "boundary", "both", "5prime"), ifelse(S2_bp_type == "boundary", "3prime", "none")), #determine if bp matches with known exon_boundaries
		bp1_frame = S1_frame, bp2_frame = S2_frame,	frame = determine_frame(S1_frame, S2_frame),
		#wt1_sequence = S1_seq, wt1_sequence_bp = S1_pos_rel, wt2_sequence = S2_seq, wt2_sequence_bp = S2_pos_rel, #determine wildtype transcripts and breakpoint positions
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
		peptide_sequence = pep_translate(cds_sequence), peptide_sequence_bp = cds_sequence_bp / 3) %>% #translate cds to peptide
	mutate(neo_peptide_sequence = ifelse(frame == "in_frame", 
		substr(peptide_sequence, round(peptide_sequence_bp - 13.5, 0), round(peptide_sequence_bp + 13.5, 0)),
		substr(peptide_sequence, round(peptide_sequence_bp - 13.5, 0), nchar(peptide_sequence))), #determine sequence potentially immunogenic
		neo_peptide_sequence_bp = peptide_sequence_bp - round(peptide_sequence_bp - 13.5, 0))
		


# select relevant columns for context_seq.csv
overlap2 = overlap %>%
	mutate(FTID = paste0(S1_symbol, "_", S1_chr, ":", S1_pos, ":", S1_strand,"_", S1_transcript_id,  "_", S2_symbol, "_", S2_chr, ":", S2_pos, ":", S2_strand, "_", S2_transcript_id)) %>%
	select(FGID, Fusion_Gene, Breakpoint1, Breakpoint2, FTID, context_sequence_id, context_sequence_100_id,
    type, exon_boundary1, exon_boundary2, exon_boundary, bp1_frame, bp2_frame, frame, 
    context_sequence, context_sequence_bp, wt1_context_sequence, wt1_context_sequence_bp, wt2_context_sequence, wt2_context_sequence_bp, 
	neo_peptide_sequence, neo_peptide_sequence_bp)
	
write.csv2(overlap2, out_file, row.names = FALSE, quote = FALSE)

# urla: write context seqs to a fasta file with the Biostrings library function writeXStringSet and write context seq based bed file for classification with R base functions
# urla: see comment in original "custom_transcriptome.py on the "problem" with bed file generation

# create vectors of context sequences and names
nameList <- c()
seqList <- c()
# create an empty dataframe for the bed file info
tmpdf <- data.frame(context_id=character(),
                    startA=numeric(),
                    endA=numeric(),
                    gene=character(),
                    emptyOnly=character(),
                    strand=character()
                    )
# each fgid+context_sequence_id needs to be processed only once, hence already processed ids are stored for reference in a list
checkList <- c()

for(i in 1:nrow(overlap2)) {
  # get names and seqs for the fasta file
#  nameList <- c(nameList, paste(overlap2$FGID[i], overlap2$context_sequence_id[i], "ft", sep = "_"))
#  nameList <- c(nameList, paste(overlap2$FGID[i], overlap2$context_sequence_id[i], "wt1", sep = "_"))
#  nameList <- c(nameList, paste(overlap2$FGID[i], overlap2$context_sequence_id[i], "wt2", sep = "_"))
  nameList <- c(nameList, paste(overlap2$FTID[i], overlap2$context_sequence_id[i], overlap2$context_sequence_bp[i], "ft", sep = "_"))
  nameList <- c(nameList, paste(overlap2$FTID[i], overlap2$context_sequence_id[i], "wt1", sep = "_"))
  nameList <- c(nameList, paste(overlap2$FTID[i], overlap2$context_sequence_id[i], "wt2", sep = "_"))
  seqList <- c(seqList, overlap2$context_sequence[i])
  seqList <- c(seqList, overlap2$wt1_context_sequence[i])
  seqList <- c(seqList, overlap2$wt2_context_sequence[i])
  
  fusion_id <- paste(overlap2$FGID[i], overlap2$context_sequence_id[i], sep = "_")
  if(! fusion_id %in% checkList) {
    # get names, length and strand for the bed file
    #ft1
    tmpdf[nrow(tmpdf) + 1, ] <- c(paste(fusion_id, "ft", sep = "_"),
                                  0, overlap2$context_sequence_bp[i], strsplit(overlap2$Fusion_Gene[i], "_")[[1]][1],
                                  "", strsplit(overlap2$Breakpoint1[i], ":")[[1]][3])
    #ft2
    tmpdf[nrow(tmpdf) + 1, ] <- c(paste(fusion_id, "ft", sep = "_"),
                                  overlap2$context_sequence_bp[i], nchar(overlap2$context_sequence[i]), strsplit(overlap2$Fusion_Gene[i], "_")[[1]][2],
                                  "", strsplit(overlap2$Breakpoint2[i], ":")[[1]][3])
    #wt11
    tmpdf[nrow(tmpdf) + 1, ] <- c(paste(fusion_id, "wt1", sep = "_"),
                                  0, overlap2$wt1_context_sequence_bp[i], strsplit(overlap2$Fusion_Gene[i], "_")[[1]][1],
                                  "", strsplit(overlap2$Breakpoint1[i], ":")[[1]][3])
    #wt12
    tmpdf[nrow(tmpdf) + 1, ] <- c(paste(fusion_id, "wt1", sep = "_"),
                                  overlap2$wt1_context_sequence_bp[i], nchar(overlap2$wt1_context_sequence[i]), strsplit(overlap2$Fusion_Gene[i], "_")[[1]][2],
                                  "", strsplit(overlap2$Breakpoint2[i], ":")[[1]][3])
    #wt21
    tmpdf[nrow(tmpdf) + 1, ] <- c(paste(fusion_id, "wt2", sep = "_"),
                                  0, overlap2$wt2_context_sequence_bp[i], strsplit(overlap2$Fusion_Gene[i], "_")[[1]][1],
                                  "", strsplit(overlap2$Breakpoint1[i], ":")[[1]][3])
    #wt22
    tmpdf[nrow(tmpdf) + 1, ] <- c(paste(fusion_id, "wt2", sep = "_"),
                                  overlap2$wt2_context_sequence_bp[i], nchar(overlap2$wt2_context_sequence[i]), strsplit(overlap2$Fusion_Gene[i], "_")[[1]][2],
                                  "", strsplit(overlap2$Breakpoint2[i], ":")[[1]][3])
    checkList <- c(checkList, fusion_id)
  }
}

# complete bed file
tmpdf$startB <- tmpdf$startA
tmpdf$endB <- tmpdf$endA
tmpdf$zeroOnly <- 0
tmpdf$oneOnly <- 1
tmpdf$startC <- -1
tmpdf[tmpdf$startA == 0, ]$startC <- tmpdf[tmpdf$startA == 0, ]$endA
tmpdf[tmpdf$startA != 0, ]$startC <- as.numeric(tmpdf[tmpdf$startA != 0, ]$endA) - as.numeric(tmpdf[tmpdf$startA != 0, ]$startA)
tmpdf$endC <- tmpdf$startA

#nameList_v2 <- paste(rep(paste(overlap2$FGID[1:nrow(overlap2)], overlap2$context_sequence_id[1:nrow(overlap2)], sep = "_"), 6), c("ft", "ft", "wt1", "wt1", "wt2", "wt2"), sep = "_")

# create an XString object with sequences, append names and write to file
seqSet <- DNAStringSet(seqList)
names(seqSet) <- nameList
writeXStringSet(seqSet, paste(out_file, ".fasta", sep = ""), append = F, format = "fasta", width = 60)
write.table(paste(length(seqSet), sum(width(seqSet)), sep = "\n"), paste(out_file, ".fasta.info", sep = ""), row.names = F, col.names = F, quote = F)
# write bed file
write.table(tmpdf, paste(out_file, ".bed", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)
print(paste0("->finished, overall processing time: ", Sys.time() - tmp_time))
