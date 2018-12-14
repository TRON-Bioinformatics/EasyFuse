#load libraries
library("optparse")

options(stringsAsFactors = FALSE)

argument_list <- list(
	make_option(c("-i", "--detected_fusions"), default="", help="Input list of detected fusions"), 
	make_option(c("-c", "--context_seq"), default="", help="Annotated context sequence list of detected fusions"), 
	make_option(c("-q", "--quantification"), default="", help="Quantification values for all fusion genes"), 
	make_option(c("-f", "--fastqc"), default="", help="FastQC file"), 
	make_option(c("-t", "--tool_state"), default="", help="tool state file"),
	make_option(c("-m", "--model"), default="", help="File containing the machine learning model"),
	make_option(c("-o", "--output"), default="", help="Final Output file for predicted fusion genes"),
)
opt <- parse_args(OptionParser(option_list=argument_list))

# check mandatory arguments
if(is.na(opt$detected_fusions) | opt$detected_fusions == "") {
	print("Mandatory parameter \"detected fusions file\" missing. Aborting...")
}
if(is.na(opt$context_seq) | opt$context_seq == "") {
	print("Mandatory parameter \"context sequence file\" missing. Aborting...")
}
if(is.na(opt$quantification) | opt$quantification == "") {
	print("Mandatory parameter \"quantification file\" missing. Aborting...")
}
if(is.na(opt$fastqc) | opt$fastqc == "") {
	print("Mandatory parameter \"fast qc file\" missing. Aborting...")
}
if(is.na(opt$tool_state) | opt$tool_state == "") {
	print("Mandatory parameter \"tool state file\" missing. Aborting...")
}
if(is.na(opt$model) | opt$model == "") {
	print("Mandatory parameter \"model file\" missing. Aborting...")
}
if(is.na(opt$output) | opt$output == "") {
	print("Mandatory parameter \"output file\" missing. Aborting...")
}

#parameter in opt structure 
opt$detected_fusions
