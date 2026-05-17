//
// EasyFuse: RANDOM_FOREST_CLASSIFIER - run the random forest classifier on detected fusions to get a set of pass fusion detections
//
//

include { PREDICTION } from '../../modules/prediction/main'

workflow RANDOM_FOREST_CLASSIFIER {

    take:
    ch_merged_results   // channel: [ val(meta), [ merged_fusions ] ]
    ch_prediction_model // channel: [prediction_model]
    ch_model_threshold  // channel: [val(threshold)]

    main:

    ch_versions = channel.empty()

    ch_prediction_input = ch_merged_results
                            .combine(ch_prediction_model)
                            .combine(ch_model_threshold)
    PREDICTION( ch_prediction_input )
    ch_versions = ch_versions.mix(PREDICTION.out.versions)

    emit:

    predictions = PREDICTION.out.predictions // channel: [ val(meta), [ fusions.pass.csv ] ]

    versions    = ch_versions                // channel: [ versions.yml ]
}
