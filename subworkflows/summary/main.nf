//
// EasyFuse: SUMMARY - summarize results from fusion prediction tools
//
//

include { MERGE_DATA } from '../../modules/utility/mergedata/main'

workflow SUMMARY {

    take:

    ch_detected_fusions    // channel: [ val(meta), [ Detected_fusions.csv ] ]
    ch_annotated_fusions   // channel: [ val(meta), annot_fusion.csv, *.debug, *.fasta ]
    ch_counts              // channel: [ val(meta), quantification.tsv ]
    ch_read_stats          // channel: [ val(meta), read_stats.tsv ]

    main:

    ch_versions = channel.empty()

    ch_mergedata_input = ch_detected_fusions
                            .join(ch_annotated_fusions)
                            .join(ch_counts)
                            .join(ch_read_stats)
    MERGE_DATA( ch_mergedata_input )

    emit:

    merged_results = MERGE_DATA.out.merged_results // channel: [val(meta), "fusion.csv"]

    versions       = ch_versions                   // channel: [ versions.yml ]
}
