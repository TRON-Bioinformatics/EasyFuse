#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC ; EASYFUSE_QC_PARSER ; EASYFUSE_SKEWER } from './modules/01_qc'
include { STAR ; EASYFUSE_READ_FILTER ; BAM2FASTQ } from './modules/02_alignment'
include { FUSION_CATCHER ; STAR_FUSION ; FUSION_CATCHER_INDEX } from './modules/03_fusion_callers'
include { EASYFUSE_FUSION_PARSER ; EASYFUSE_FUSION_ANNOTATION } from './modules/04_joint_fusion_calling'
include { EASYFUSE_REQUANTIFY_STAR_INDEX ; EASYFUSE_REQUANTIFY_FILTER ; 
    EASYFUSE_REQUANTIFY_BAM2FASTQ ; EASYFUSE_REQUANTIFY_STAR ; 
    EASYFUSE_REQUANTIFY_COUNT } from './modules/05_requantification'
include { EASYFUSE_SUMMARIZE_DATA } from './modules/06_summarize'
	    

def helpMessage() {
    log.info params.help_message
}

if (params.help) {
    helpMessage()
    exit 0
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

workflow TOOLS {
    take:
    filtered_fastq
    chimeric_reads

    main:
    FUSION_CATCHER(ALIGNMENT.out.fastqs)
    STAR_FUSION(ALIGNMENT.out.chimeric_reads)

    emit:
    fusioncatcher_results = FUSIONCATCHER.out.fusions
    starfusion_results = STAR_FUSION.out.fusions
}

workflow ANNOTATION {
    take:
    fusioncatcher_result
    starfusion_result

    main:
    EASYFUSE_FUSION_PARSER(
        fusioncatcher_result.join(
        starfusion_result
    ))
    EASYFUSE_FUSION_ANNOTATION(EASYFUSE_FUSION_PARSER.out.fusions)

    emit:
    annotated_fusions = EASYFUSE_FUSION_ANNOTATION.out.annot_fusions
}

workflow REQUANT {
    take:
    filtered_fastq

    main:
    EASYFUSE_REQUANTIFY_STAR_INDEX(EASYFUSE_FUSION_ANNOTATION.out.annot_fusions)
    EASYFUSE_REQUANTIFY_FILTER(
        ALIGNMENT.out.bams.join(
        EASYFUSE_FUSION_ANNOTATION.out.annot_fusions).join(
        ALIGNMENT.out.read_stats
    ))
    EASYFUSE_REQUANTIFY_BAM2FASTQ(EASYFUSE_REQUANTIFY_FILTER.out.bams)
    EASYFUSE_REQUANTIFY_STAR(
        EASYFUSE_REQUANTIFY_BAM2FASTQ.out.fastqs.join(
        EASYFUSE_REQUANTIFY_STAR_INDEX.out.star_index
    ))
    EASYFUSE_REQUANTIFY_COUNT(
        EASYFUSE_REQUANTIFY_STAR.out.bams.join(
        EASYFUSE_REQUANTIFY_STAR.out.read_stats
    ))

    emit:
    requant_results = EASYFUSE_REQUANTIFY_COUNT.out.counts
}


workflow {

    QC(input_files)
    ALIGNMENT(QC.out.trimmed_fastq)

    // fusion calling
    FUSION_CATCHER(ALIGNMENT.out.fastqs)
    STAR_FUSION(ALIGNMENT.out.chimeric_reads)

    // fusion merging
    EASYFUSE_FUSION_PARSER(
        FUSION_CATCHER.out.fusions.join(
        STAR_FUSION.out.fusions
    ))
    EASYFUSE_FUSION_ANNOTATION(EASYFUSE_FUSION_PARSER.out.fusions)

    // requantification
    EASYFUSE_REQUANTIFY_STAR_INDEX(EASYFUSE_FUSION_ANNOTATION.out.annot_fusions)
    EASYFUSE_REQUANTIFY_FILTER(
        ALIGNMENT.out.bams.join(
        EASYFUSE_FUSION_ANNOTATION.out.annot_fusions).join(
        ALIGNMENT.out.read_stats
    ))
    EASYFUSE_REQUANTIFY_BAM2FASTQ(EASYFUSE_REQUANTIFY_FILTER.out.bams)
    EASYFUSE_REQUANTIFY_STAR(
        EASYFUSE_REQUANTIFY_BAM2FASTQ.out.fastqs.join(
        EASYFUSE_REQUANTIFY_STAR_INDEX.out.star_index
    ))
    EASYFUSE_REQUANTIFY_COUNT(
        EASYFUSE_REQUANTIFY_STAR.out.bams.join(
        EASYFUSE_REQUANTIFY_STAR.out.read_stats
    ))

    EASYFUSE_SUMMARIZE_DATA(
        EASYFUSE_FUSION_PARSER.out.fusions.join(
        EASYFUSE_FUSION_ANNOTATION.out.annot_fusions).join(
        EASYFUSE_REQUANTIFY_COUNT.out.counts).join(
        EASYFUSE_REQUANTIFY_STAR.out.read_stats
    ))
}
