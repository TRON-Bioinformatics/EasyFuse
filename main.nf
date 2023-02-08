#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC ; EASYFUSE_QC_PARSER ; EASY_FUSE_SKEWER } from './modules/01_qc'

params.help= false
params.input_files = false

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
    EASY_FUSE_SKEWER(input_files.join(EASYFUSE_QC_PARSER.out.qc_table))

}
