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

    main:

    ch_versions = channel.empty()


    ch_fusionparser_input = channel.empty()
                                .mix(ch_fusioncatcher_results.ifEmpty { [[]] })
                                .mix(ch_starfusion_results.ifEmpty { [[]] })
                                .mix(ch_arriba_results.ifEmpty { [[]] })
                                .groupTuple()
                                .map { meta, fusioncatcher_fusions, starfusion_fusions, arriba_fusions ->
                                    tuple(meta, fusioncatcher_fusions, starfusion_fusions, arriba_fusions)
                                }
    FUSION_PARSER ( ch_fusionparser_input )
    ch_versions = ch_versions.mix(FUSION_PARSER.out.versions)

    FUSIONANNOTATER ( FUSION_PARSER.out.fusions )
    ch_versions = ch_versions.mix(FUSIONANNOTATER.out.versions)

    emit:

    annot_fusions    = FUSIONANNOTATER.out.annot_fusions // channel: [ val(meta), [${prefix}_Annotated_fusions.csv]]
    detected_fusions = FUSION_PARSER.out.fusions         // channel: [ val(meta), [${prefix}_Detected_fusions.csv]]

    versions = ch_versions                              // channel: [ versions.yml ]
}
