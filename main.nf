#!/usr/bin/env nextflow

include { EASYFUSE         } from './workflows/easyfuse'
include { INPUT_VALIDATION } from './subworkflows/validation/parameter_validation'

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow TRONBIOINFORMATICS_EASYFUSE {

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
    versions = EASYFUSE.out.versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    INPUT_VALIDATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.fusion_tools,
        params.model_pred,
        params.model_prefix,
        params.model_threshold,
        params.reference,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    TRONBIOINFORMATICS_EASYFUSE (
        INPUT_VALIDATION.out.samplesheet,
        INPUT_VALIDATION.out.fusiontools,
        INPUT_VALIDATION.out.reference_fasta,
        INPUT_VALIDATION.out.reference_gtf,
        INPUT_VALIDATION.out.reference_tsl,
        INPUT_VALIDATION.out.annotation_db,
        INPUT_VALIDATION.out.starfusion_index,
        INPUT_VALIDATION.out.fusioncatcher_index,
        INPUT_VALIDATION.out.stararriba_index,
        INPUT_VALIDATION.out.prediction_model,
        INPUT_VALIDATION.out.model_threshold
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
