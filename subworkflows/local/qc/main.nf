//
// EasyFuse: QC - perform pre-alignment qc and trimming on raw rna-seq reads
//
//

include { FASTP } from '../../../modules/local/fastp/main'

workflow QC {

    take:
    ch_fastqs // channel: [ val(meta), fastq1, fastq2 ]

    main:

    ch_versions = channel.empty()

    FASTP ( ch_fastqs )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    emit:

    trimmed_fastqs = FASTP.out.trimmed_fastqs  // channel: [ val(meta), trimmed_fastq1, trimmed_fastq2 ]
    versions       = ch_versions               // channel: [ versions.yml ]
}
