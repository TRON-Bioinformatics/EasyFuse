#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC ; EASYFUSE_QC_PARSER ; EASYFUSE_SKEWER } from './modules/01_qc'
include { STAR ; EASYFUSE_READ_FILTER ; BAM2FASTQ } from './modules/02_alignment'
include { MAPSPLICE ; FUSION_CATCHER ; STAR_FUSION} from './modules/03_fusion_callers'

params.help= false
params.input_files = false
params.reference = false
params.chromosome_dir = "/projects/data/human/ensembl/GRCh38.86/fasta"
params.bowtie_index = "/projects/data/human/ensembl/GRCh38.86/bowtie_index/hg38"
params.gtf = "/projects/data/human/ensembl/GRCh38.86/Homo_sapiens.GRCh38.86.gtf"
params.fusioncatcher_index = "/scratch/info/data/easyfuse/easyfuse_ref/fusioncatcher_index/"
params.starfusion_index = "/projects/data/human/ensembl/GRCh38.86/starfusion_index/"
params.output = false

def helpMessage() {
    log.info params.help_message
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.reference) {
    log.error "--reference is required"
    exit 1
}

if (!params.output) {
    log.error "--output is required"
    exit 1
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'fastq1', 'fastq2'], sep: "\t")
    .map{ row-> tuple(row.name, row.fastq1, row.fastq2) }
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}


workflow {

    // QC
    FASTQC(input_files)
    EASYFUSE_QC_PARSER(FASTQC.out.qc_data)
    EASYFUSE_SKEWER(input_files.join(EASYFUSE_QC_PARSER.out.qc_table))

    // Alignment
    STAR(EASYFUSE_SKEWER.out.trimmed_fastq)
    EASYFUSE_READ_FILTER(STAR.out.bams)
    BAM2FASTQ(EASYFUSE_READ_FILTER.out.bams)

    // Fusions
    MAPSPLICE(BAM2FASTQ.out.fastqs)
    FUSION_CATCHER(BAM2FASTQ.out.fastqs)
    STAR_FUSION(STAR.out.chimeric_reads)
}
