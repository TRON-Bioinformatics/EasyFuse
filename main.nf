#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC ; EASYFUSE_QC_PARSER ; EASYFUSE_SKEWER } from './modules/01_qc'
include { STAR ; EASYFUSE_READ_FILTER ; BAM2FASTQ } from './modules/02_alignment'
include { FUSION_CATCHER ; STAR_FUSION ; FUSION_CATCHER_INDEX } from './modules/03_fusion_callers'
include { EASYFUSE_FUSION_PARSER ; EASYFUSE_FUSION_ANNOTATION ; EASYFUSE_STAR_INDEX ;
            EASYFUSE_REQUANTIFY_FILTER} from './modules/04_joint_fusion_calling'

params.help = false
params.index = false
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

if (!params.index && !params.reference) {
    log.error "--reference is required"
    exit 1
}

if (!params.index && !params.output) {
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
} else if (!params.index) {
  exit 1, "Input file not specified!"
}

workflow QC {
    take:
    input_files

    main:
    FASTQC(input_files)
    EASYFUSE_QC_PARSER(FASTQC.out.qc_data)
    EASYFUSE_SKEWER(input_files.join(EASYFUSE_QC_PARSER.out.qc_table))

    emit:
    trimmed_fastq = EASYFUSE_SKEWER.out.trimmed_fastq
}

workflow ALIGNMENT {
    take:
    trimmed_fastq

    main:
    STAR(trimmed_fastq)
    EASYFUSE_READ_FILTER(STAR.out.bams)
    BAM2FASTQ(EASYFUSE_READ_FILTER.out.bams)

    emit:
    chimeric_reads = STAR.out.chimeric_reads
    fastqs = BAM2FASTQ.out.fastqs
    bams = EASYFUSE_READ_FILTER.out.bams
    read_stats = STAR.out.read_stats
}


workflow {

    if (params.index) {
        FUSION_CATCHER_INDEX()
    } else {
        QC(input_files)

        ALIGNMENT(QC.out.trimmed_fastq)

        // Fusions
        //FUSION_CATCHER(ALIGNMENT.out.fastqs)
        STAR_FUSION(ALIGNMENT.out.chimeric_reads)

        // joint fusion calling
        //EASYFUSE_FUSION_PARSER(FUSION_CATCHER.out.fusions.join(STAR_FUSION.out.fusions))
        EASYFUSE_FUSION_PARSER(STAR_FUSION.out.fusions)
        EASYFUSE_FUSION_ANNOTATION(EASYFUSE_FUSION_PARSER.out.fusions)
        EASYFUSE_STAR_INDEX(EASYFUSE_FUSION_ANNOTATION.out.fusions)
        EASYFUSE_REQUANTIFY_FILTER(ALIGNMENT.out.bams.join(
            EASYFUSE_FUSION_ANNOTATION.out.fusions_long).join(
            ALIGNMENT.out.read_stats))
    }
}
