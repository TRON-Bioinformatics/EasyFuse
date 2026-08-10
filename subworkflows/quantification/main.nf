//
// EasyFuse: REQUANTIFICATION - subworkflow for requantifying fusion events using BPQuant
//
//

include { FUSION_FILTER     } from '../../modules/utility/fusionfilter/main'
include { BAM2FASTQ         } from '../../modules/bam2fastq/main'
include { FUSION2CSV        } from '../../modules/utility/fusion2csv/main'
include { BPQUANT_CSV2FASTA } from '../../modules/bpquant/csv2fasta/main'
include { BPQUANT_INDEX     } from '../../modules/bpquant/index/main'
include { BPQUANT_ALIGN     } from '../../modules/bpquant/align/main'
include { BPQUANT_COUNT     } from '../../modules/bpquant/count/main'


workflow QUANTIFICATION {

    take:

    ch_bam                // channel: [ val(meta), [ bam ] ]
    ch_read_stats         // channel: [ val(meta), [read_stats] ]
    ch_annotated_fusions  // channel: [ val(meta), [annotated_fusions] ]

    main:

    ch_versions = channel.empty()

    ch_fusionsfilter_input = ch_bam
                                .join(ch_annotated_fusions)
                                .join(ch_read_stats)
                                .map { meta, bam, annot_fusions_csv, annot_fusions_debug, annot_fusions_fasta, read_stats
                                    -> tuple(meta, bam, annot_fusions_csv, annot_fusions_debug, annot_fusions_fasta, read_stats)
                                }
    FUSION_FILTER ( ch_fusionsfilter_input )

    BAM2FASTQ ( FUSION_FILTER.out.bam )
    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions)

    FUSION2CSV (ch_annotated_fusions)

    BPQUANT_CSV2FASTA (FUSION2CSV.out.formatted_csv)
    ch_versions = ch_versions.mix(BPQUANT_CSV2FASTA.out.versions)

    BPQUANT_INDEX (BPQUANT_CSV2FASTA.out.formatted_fasta)
    ch_versions = ch_versions.mix(BPQUANT_INDEX.out.versions)

    BPQUANT_ALIGN (BAM2FASTQ.out.fastqs.combine(BPQUANT_INDEX.out.star_index, by: 0))
    ch_versions = ch_versions.mix(BPQUANT_ALIGN.out.versions)

    BPQUANT_COUNT (BPQUANT_ALIGN.out.bam.join(FUSION2CSV.out.formatted_csv))
    ch_versions = ch_versions.mix(BPQUANT_COUNT.out.versions)

    emit:

    counts     = BPQUANT_COUNT.out.counts      // channel: [ val(meta), path(quantification.tsv) ]
    read_stats = BPQUANT_ALIGN.out.read_stats  // channel: [ val(meta), path(read_stats.tsv) ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
