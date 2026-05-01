//
// EasyFuse: Alignment - align paired end RNA-seq reads to a reference genome and transcriptome,
// and filter for fusion-supporting reads
//

include { STAR_ALIGN  } from '../../../modules/local/star/main'
include { BAM2FASTQ   } from '../../../modules/local/bam2fastq/main'
include { READ_FILTER } from '../../../modules/local/utility/readfilter/main'

workflow READ_FILTERING {

    take:
    ch_trimmed_fastqs // channel: [ val(meta), [ trimmed_fastq1, trimmed_fastq2 ] ]

    main:

    ch_versions = channel.empty()

    STAR_ALIGN ( ch_trimmed_fastqs )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    READ_FILTER ( STAR_ALIGN.out.bam )
    ch_versions = ch_versions.mix(READ_FILTER.out.versions)

    BAM2FASTQ ( READ_FILTER.out.filtered_bam )
    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions)

    emit:
    bam            = READ_FILTER.out.filtered_bam    // channel: [ val(meta), [ bam ] ]
    fastqs         = BAM2FASTQ.out.fastqs            // channel: [ val(meta), [ fastq1, fastq2 ] ]
    read_stats     = STAR_ALIGN.out.read_stats       // channel: [ val(meta), [ read_stats ] ]
    chimeric_reads = STAR_ALIGN.out.chimeric_reads   // channel: [ val(meta), [ chimeric_reads ] ]

    versions       = ch_versions                     // channel: [ versions.yml ]
}
