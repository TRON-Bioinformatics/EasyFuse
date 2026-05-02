//
// EasyFuse: FUSION_ANNOTATION - subworkflow to annotate gene fusions detected using tools of choice
//
//

include { FUSION_PARSER   } from '../../../modules/local/fusionparsing/fusionparser/main'
include { FUSIONANNOTATER } from '../../../modules/local/fusionannotation/main'

workflow FUSION_ANNOTATION {

    take:
    ch_fusioncatcher_results // channel: [ val(meta), [ fusioncatcher_fusions ] ]
    ch_starfusion_results    // channel: [ val(meta), [ starfusion_fusions ]]
    ch_arriba_results        // channel: [ val(meta), [ arriba_fusions ]]
    ch_annotation_db         // channel: annotation db
    ch_reference_fasta       // channel: reference fasta
    ch_reference_tsl         // channel: reference tsl
    ch_fusiontools           // [ [arriba: true, fusioncatcher: true, starfusion: true] ]

    main:

    ch_versions = channel.empty()

    // parse detected fusions
    ch_arriba = ch_fusiontools.value.run_arriba ? ch_arriba_results : channel.empty()
    ch_starfusion = ch_fusiontools.value.run_starfusion ? ch_starfusion_results : channel.empty()
    ch_fusioncatcher = ch_fusiontools.value.run_fusioncatcher ? ch_fusioncatcher_results : channel.empty()

    // generate fusionparser input, based on the tools run
    ch_fusionparser_input = channel
                                .empty()
                                .mix(ch_fusioncatcher.map { meta, value -> tuple(meta, value, null, null) })
                                .mix(ch_starfusion.map { meta, value -> tuple(meta, null, value, null) })
                                .mix(ch_arriba.map { meta, value -> tuple(meta, null, null, value) })
                                .groupTuple()
                                .map { meta, fc_list, sf_list, ar_list ->

                                    def arriba = ar_list.find { it -> it != null }
                                    def starfusion = sf_list.find { it -> it != null }
                                    def fusioncatcher = fc_list.find { it -> it != null }

                                    tuple(
                                        meta,
                                        fusioncatcher ?: [],
                                        starfusion ?: [],
                                        arriba ?: []
                                    )
                                }
    FUSION_PARSER ( ch_fusionparser_input )
    ch_versions = ch_versions.mix(FUSION_PARSER.out.versions)

    ch_fusionannotater_input = FUSION_PARSER.out.fusions
                                .combine(ch_annotation_db)
                                .combine(ch_reference_fasta)
                                .combine(ch_reference_tsl)
                                .map { meta, fusions, annotation_db, ref_fasta, ref_tsl ->
                                    tuple(meta, fusions, annotation_db, ref_fasta, ref_tsl)
                                }
    FUSIONANNOTATER ( ch_fusionannotater_input )
    ch_versions = ch_versions.mix(FUSIONANNOTATER.out.versions)

    emit:

    annot_fusions    = FUSIONANNOTATER.out.annot_fusions // channel: [ val(meta), [${prefix}_Annotated_fusions.csv]]
    detected_fusions = FUSION_PARSER.out.fusions         // channel: [ val(meta), [${prefix}_Detected_fusions.csv]]

    versions = ch_versions                               // channel: [ versions.yml ]
}
