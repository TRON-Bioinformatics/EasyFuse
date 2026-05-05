#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRON/easyfuse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/TRON/easyfuse
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EASYFUSE                } from './workflows/easyfuse'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_easyfuse_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_easyfuse_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_easyfuse_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow TRONPRIVATE_EASYFUSE {

    take:
    samplesheet             // channel: samplesheet read in from --input
    fusion_tools            // channel: fusion tools to run (read in from params.fusion_tools)
    ch_reference_fasta      // channel: [reference fasta (read in from --reference)]
    ch_reference_gtf        // channel: [reference gtf (read in from --reference)]
    ch_reference_tsl        // channel: [reference tsl (read in from --reference)]
    ch_annotation_db        // channel: [annotation db]
    ch_starfusion_index     // channel: [starfusion index]
    ch_fusioncatcher_index  // channel: [fusioncatcher index]
    ch_stararriba_index     // channel: [stararriba index]
    ch_prediction_model     // channel: [prediction model]
    ch_model_threshold      // channel: [val(threshold)]

    main:

    //
    // WORKFLOW: Run pipeline
    //
    EASYFUSE (
        samplesheet,
        fusion_tools,
        ch_reference_fasta,
        ch_reference_gtf,
        ch_reference_tsl,
        ch_annotation_db,
        ch_starfusion_index,
        ch_fusioncatcher_index,
        ch_stararriba_index,
        ch_prediction_model,
        ch_model_threshold
    )
    emit:
    multiqc_report = EASYFUSE.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.fusion_tools,
        params.ensembl_version,
        params.model_pred,
        params.model_threshold,
        params.reference,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    TRONPRIVATE_EASYFUSE (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.fusiontools,
        PIPELINE_INITIALISATION.out.reference_fasta,
        PIPELINE_INITIALISATION.out.reference_gtf,
        PIPELINE_INITIALISATION.out.reference_tsl,
        PIPELINE_INITIALISATION.out.annotation_db,
        PIPELINE_INITIALISATION.out.starfusion_index,
        PIPELINE_INITIALISATION.out.fusioncatcher_index,
        PIPELINE_INITIALISATION.out.stararriba_index,
        PIPELINE_INITIALISATION.out.prediction_model,
        PIPELINE_INITIALISATION.out.model_threshold
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
        TRONPRIVATE_EASYFUSE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
