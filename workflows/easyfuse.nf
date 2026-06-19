/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/validation/pipeline_utils'
include { softwareVersionsToYAML } from '../subworkflows/validation/pipeline_utils'

include { QC                        } from '../subworkflows/qc/main'
include { READ_FILTERING            } from '../subworkflows/read_filtering/main'
include { FUSION_PREDICTION         } from '../subworkflows/fusion_prediction/main'
include { FUSION_ANNOTATION         } from '../subworkflows/fusion_annotation/main'
include { QUANTIFICATION            } from '../subworkflows/quantification/main'
include { SUMMARY                   } from '../subworkflows/summary/main'
include { RANDOM_FOREST_CLASSIFIER  } from '../subworkflows/random_forest_classifier/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EASYFUSE {

    take:

    ch_samplesheet             // channel: samplesheet read in from --input
    ch_fusiontools             // channel: [ val([run_arriba: true, run_fusioncatcher: true, run_starfusion: true]) ]
    ch_reference_fasta         // channel: reference fasta (read in from --reference)
    ch_reference_gtf           // channel: reference gtf (read in from --reference)
    ch_reference_tsl           // channel: reference tsl (read in from --reference)
    ch_annotation_db           // channel: annotation db
    ch_starfusion_index        // channel: starfusion index
    ch_fusioncatcher_index     // channel: fusioncatcher index
    ch_stararriba_index        // channel: stararriba index
    ch_prediction_model        // channel: [prediction model]
    ch_model_threshold         // channel: [val(threshold)]

    main:

    ch_versions            = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        QC Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    QC ( ch_samplesheet )
    ch_versions = ch_versions.mix(QC.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Read-Filtering Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    READ_FILTERING (
        QC.out.trimmed_fastqs,
        ch_stararriba_index
    )
    ch_versions = ch_versions.mix(READ_FILTERING.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Fusion-Prediction Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    FUSION_PREDICTION (
        READ_FILTERING.out.fastqs,
        ch_fusioncatcher_index,
        ch_starfusion_index,
        ch_stararriba_index,
        ch_reference_gtf,
        ch_reference_fasta,
        ch_fusiontools
    )
    ch_versions = ch_versions.mix(FUSION_PREDICTION.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Fusion-Annotation Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    FUSION_ANNOTATION (
        FUSION_PREDICTION.out.fusioncatcher_results,
        FUSION_PREDICTION.out.starfusion_results,
        FUSION_PREDICTION.out.arriba_results,
        ch_annotation_db,
        ch_reference_fasta,
        ch_reference_tsl,
        ch_fusiontools
    )
    ch_versions = ch_versions.mix(FUSION_ANNOTATION.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Quantification Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    QUANTIFICATION (
        READ_FILTERING.out.bam,
        READ_FILTERING.out.read_stats,
        FUSION_ANNOTATION.out.annot_fusions
    )
    ch_versions = ch_versions.mix(QUANTIFICATION.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Summary Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    SUMMARY (
        FUSION_ANNOTATION.out.detected_fusions,
        FUSION_ANNOTATION.out.annot_fusions,
        QUANTIFICATION.out.counts,
        QUANTIFICATION.out.read_stats,
        ch_fusiontools
    )
    ch_versions = ch_versions.mix(SUMMARY.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Random Forest Classifier Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    RANDOM_FOREST_CLASSIFIER (
        SUMMARY.out.merged_results,
        ch_prediction_model,
        ch_model_threshold
    )
    ch_versions = ch_versions.mix(RANDOM_FOREST_CLASSIFIER.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collate and save software versions
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/version_info",
            name:  'easyfuse_software_' + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:

    fusions        = SUMMARY.out.merged_results                // channel: [path(fusions.csv)]
    fusions_pass   = RANDOM_FOREST_CLASSIFIER.out.predictions  // channel: [path(fusions.pass.csv)]

    versions       = ch_versions                               // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
