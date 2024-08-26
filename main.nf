#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTP } from './modules/01_qc'
include { STAR ; STAR_ARRIBA ; READ_FILTER ; BAM2FASTQ } from './modules/02_alignment'
include { FUSION_CATCHER ; STAR_FUSION ; ARRIBA } from './modules/03_fusion_callers'
include { FUSION_PARSER ; FUSION_ANNOTATION } from './modules/04_joint_fusion_calling'
include { FUSION2CSV ; CSV2FASTA ; STAR_INDEX ; FUSION_FILTER ; STAR_CUSTOM ; READ_COUNT } from './modules/05_requantification'
include { MERGE_DATA ; PREDICTION } from './modules/06_summarize'
	    

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

if (!params.reference) {
    log.error "--reference is required"
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
    FASTP(input_files)

    emit:
    trimmed_fastq = FASTP.out.trimmed_fastq
}

workflow ALIGNMENT {
    take:
    trimmed_fastq

    main:
    STAR(trimmed_fastq)
    READ_FILTER(STAR.out.bams)
    BAM2FASTQ(READ_FILTER.out.bams)

    emit:
    chimeric_reads = STAR.out.chimeric_reads
    fastqs = BAM2FASTQ.out.fastqs
    bams = READ_FILTER.out.bams
    read_stats = STAR.out.read_stats
}

workflow TOOLS {
    take:
    filtered_fastqs
    chimeric_reads

    main:
    FUSION_CATCHER(filtered_fastqs)
    STAR_FUSION(filtered_fastqs)
    STAR_ARRIBA(filtered_fastqs)
    ARRIBA(STAR_ARRIBA.out.bams)


    emit:
    fusioncatcher_results = FUSION_CATCHER.out.fusions
    starfusion_results = STAR_FUSION.out.fusions
    arriba_results = ARRIBA.out.fusions
}

workflow ANNOTATION {
    take:
    fusioncatcher_results
    starfusion_results
    arriba_results

    main:
    FUSION_PARSER(
        fusioncatcher_results.join(
        starfusion_results.join(
        arriba_results)
    ))
    FUSION_ANNOTATION(FUSION_PARSER.out.fusions)

    emit:
    fusions = FUSION_PARSER.out.fusions
    annotated_fusions = FUSION_ANNOTATION.out.annot_fusions
}

workflow REQUANTIFICATION {
    take:
    annotated_fusions
    bams
    read_stats

    main:
    FUSION_FILTER(
        bams.join(
        annotated_fusions).join(
        read_stats
    ))
    BAM2FASTQ(FUSION_FILTER.out.bams)
    FUSION2CSV(annotated_fusions)
    CSV2FASTA(FUSION2CSV.out.formatted_csv)
    STAR_INDEX(CSV2FASTA.out.formatted_fasta)
    STAR_CUSTOM(
        BAM2FASTQ.out.fastqs.join(
        STAR_INDEX.out.star_index
    ))
    READ_COUNT(
        STAR_CUSTOM.out.bams.join(
        FUSION2CSV.out.formatted_csv
    ))

    emit:
    counts = READ_COUNT.out.counts
    read_stats = STAR_CUSTOM.out.read_stats
}


workflow {

    // quality trimming
    QC(input_files)

    // read alignment and read filtering, required as tool input
    ALIGNMENT(QC.out.trimmed_fastq)

    // fusion calling
    TOOLS(
        ALIGNMENT.out.fastqs,
        ALIGNMENT.out.chimeric_reads
    )

    // fusion merging
    ANNOTATION(
        TOOLS.out.fusioncatcher_results,
        TOOLS.out.starfusion_results,
        TOOLS.out.arriba_results
    )

    // requantification
    REQUANTIFICATION(
        ANNOTATION.out.annotated_fusions,
        ALIGNMENT.out.bams,
        ALIGNMENT.out.read_stats
    )

    // summarize results
    MERGE_DATA(
        ANNOTATION.out.fusions.join(
        ANNOTATION.out.annotated_fusions).join(
        REQUANTIFICATION.out.counts).join(
        REQUANTIFICATION.out.read_stats
    ))

    // run prediction
    PREDICTION(
        MERGE_DATA.out.merged_results
    )
}
