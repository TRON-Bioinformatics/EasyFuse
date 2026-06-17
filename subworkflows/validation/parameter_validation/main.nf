//
// Subworkflow with functionality specific to the TRON-Bioinformatics/easyfuse pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../schema_validation'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionSummary         } from '../pipeline_utils'
include { UTILS_NFCORE_PIPELINE     } from '../pipeline_utils'
include { UTILS_NEXTFLOW_PIPELINE   } from '../pipeline_utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INPUT_VALIDATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args // array: List of positional nextflow CLI args
    outdir            // string: The output directory where the results will be saved
    input             // string: Path to input samplesheet
    fusion_tools      // string: comma separated string of fusion prediction tools
    ensembl_version   // string: ensembl version info.
    model_pred        // string: path to the random forest classifier model
    model_threshold   // number: model threshold for the random forest classifier
    reference         // string: Path to reference directory containing genome files (fasta, gtf, star indices etc.)
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    channel
    .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
    .map {
        meta, fastq_1, fastq_2 ->
            if (!fastq_2) {
                return [ meta.id, meta + [ paired_end:false ],  [ fastq_1 ] ]
            } else {
                return [ meta.id, meta + [ paired_end:true ], [ fastq_1, fastq_2 ] ]
            }
    }
    .groupTuple()
    .map { samplesheet ->
        def (_id, metas, fastqs) = samplesheet
        workflow.profile.contains('test') ? [ metas[0], fastqs ] : validateInputSamplesheet(samplesheet)
    }
    .map {
        meta, fastqs ->
            return [ meta, fastqs[0], fastqs[1] ]
    }
    .set { ch_samplesheet }


    //
    // Validate fusion tools provided by the user
    //
    fusiontools = validateFusionTools(fusion_tools)

    //
    // Build reference file channels
    //

    // get reference fasta
    ch_reference_fasta = channel.value(
        file(reference.toString().replaceFirst(/\/$/, '') + '/Homo_sapiens.GRCh38.dna.primary_assembly.fa', checkIfExists: true)
    )

    // get reference gtf
    ch_reference_gtf = channel.value(
        file(reference.toString().replaceFirst(/\/$/, '') + "/Homo_sapiens.GRCh38.${ensembl_version}.gtf", checkIfExists: true)
    )

    // ch reference tsl
    ch_reference_tsl = channel.value(
        file(reference.toString().replaceFirst(/\/$/, '') + "/Homo_sapiens.GRCh38.${ensembl_version}.gtf.tsl", checkIfExists: true)
    )

    // get annotation db
    ch_annotation_db = channel.value(
        file(reference.toString().replaceFirst(/\/$/, '') + "/Homo_sapiens.GRCh38.${ensembl_version}.gff3.db", checkIfExists: true)
    )

    // get starfusion index
    ch_starfusion_index = channel.value(
        file(reference.toString().replaceFirst(/\/$/, '') + "/starfusion_index", checkIfExists: true)
    )

    // get fusioncatcher index
    ch_fusioncatcher_index = channel.value(
        file(reference.toString().replaceFirst(/\/$/, '') + "/fusioncatcher_index", checkIfExists: true)
    )

    // get stararriba index
    ch_stararriba_index = channel.value(
        file(reference.toString().replaceFirst(/\/$/, '') + "/star_index", checkIfExists: true)
    )

    // random forest classifier model
    ch_prediction_model = channel.value(
        file("${projectDir}/assets/data/model/${model_pred}", checkIfExists: true)
    )

    // model threshold
    ch_model_threshold = channel.value(model_threshold)


    emit:

    samplesheet             = ch_samplesheet
    fusiontools             = fusiontools
    reference_fasta         = ch_reference_fasta
    reference_gtf           = ch_reference_gtf
    reference_tsl           = ch_reference_tsl
    annotation_db           = ch_annotation_db
    starfusion_index        = ch_starfusion_index
    fusioncatcher_index     = ch_fusioncatcher_index
    stararriba_index        = ch_stararriba_index
    prediction_model        = ch_prediction_model
    model_threshold         = ch_model_threshold

    versions            = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
}


//
// Validate if the user provided fusion tools are a valid choice
//
def validateFusionTools(fusion_tools) {
    if (!fusion_tools) {
        error "Please provide --fusion_tools. Valid options: arriba, starfusion, fusioncatcher"
    }

    def tools = fusion_tools.split(',').collect { tool -> tool.trim().toLowerCase() }
    def valid_tools = ['arriba', 'starfusion', 'fusioncatcher']

    def invalid = tools - valid_tools
    if (invalid) {
        error "Invalid fusion tool(s): ${invalid.join(', ')}. Valid options: ${valid_tools.join(', ')}"
    }

    return [
        run_arriba       : 'arriba' in tools,
        run_starfusion   : 'starfusion' in tools,
        run_fusioncatcher: 'fusioncatcher' in tools
    ]
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs.flatten() ]  // flatten here, not in the downstream map
}
