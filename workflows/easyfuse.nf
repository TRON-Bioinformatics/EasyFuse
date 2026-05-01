/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_easyfuse_pipeline'

include { FASTP as QC                            } from '../modules/local/fastp/main'
include { READ_FILTERING                         } from '../subworkflows/local/read_filtering/main'
include { FUSION_PREDICTION                      } from '../subworkflows/local/fusion_prediction/main'
include { FUSION_ANNOTATION                      } from '../subworkflows/local/fusion_annotation/main'
include { QUANTIFICATION                         } from '../subworkflows/local/quantification/main'
include { MERGE_DATA as SUMMARY                  } from '../modules/local/utility/mergedata/main'
include { PREDICTION as RANDOM_FOREST_CLASSIFIER } from '../modules/local/prediction/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EASYFUSE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fusiontools // channel: fusion tools to run (read in from params.fusion_tools)

    main:

    ch_versions            = channel.empty()
    ch_multiqc_files       = channel.empty()
    ch_reference_gtf       = channel.empty()
    ch_reference_fasta     = channel.empty()
    ch_starfusion_index    = channel.empty()
    ch_stararriba_index    = channel.empty()
    ch_fusioncatcher_index = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        QC Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    QC( ch_samplesheet )
    ch_versions = ch_versions.mix(QC.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Read-Filtering Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    READ_FILTERING( QC.out.trimmed_fastq )
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
    ch_arriba_fusions        = FUSION_PREDICTION.out.arriba_results
    ch_starfusion_fusions    = FUSION_PREDICTION.out.starfusion_results
    ch_fusioncatcher_fusions = FUSION_PREDICTION.out.fusioncatcher_results

    FUSION_ANNOTATION (
        ch_fusioncatcher_fusions,
        ch_starfusion_fusions,
        ch_arriba_fusions
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
        FUSION_ANNOTATION.out.detected_fusions
            .join(FUSION_ANNOTATION.out.annot_fusions)
            .join(QUANTIFICATION.out.counts)
            .join(QUANTIFICATION.out.read_stats)
    )
    ch_versions = ch_versions.mix(SUMMARY.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Random Forest Classifier Layer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    RANDOM_FOREST_CLASSIFIER (
        SUMMARY.out.merged_data
    )
    ch_versions = ch_versions.mix(RANDOM_FOREST_CLASSIFIER.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collate and save software versions
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'easyfuse_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    ch_multiqc_config        = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()
    summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary      = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files         = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_files         = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions            = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
