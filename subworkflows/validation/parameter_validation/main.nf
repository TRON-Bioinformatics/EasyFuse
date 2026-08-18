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

include { UNTAR as UNTAR_STAR_INDEX          } from '../../../modules/utility/untar/main'
include { UNTAR as UNTAR_STARFUSION_INDEX    } from '../../../modules/utility/untar/main'
include { UNTAR as UNTAR_FUSIONCATCHER_INDEX } from '../../../modules/utility/untar/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INPUT_VALIDATION {

    take:
    version                 // boolean: Display version and exit
    validate_params         // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs        // boolean: Do not use coloured log outputs
    nextflow_cli_args       // array: List of positional nextflow CLI args
    _outdir                 // string: The output directory where the results will be saved
    input                   // string: Path to input samplesheet
    fusion_tools            // string: comma separated string of fusion prediction tools
    model_pred              // string: path to the random forest classifier model
    model_prefix            // string: prefix for the random forest classifier model
    model_threshold         // number: model threshold for the random forest classifier
    reference               // string: Path to reference directory containing genome files (fasta, gtf, star indices etc.)
    help                    // boolean: Display help message and exit
    help_full               // boolean: Show the full help message
    show_hidden             // boolean: Show hidden parameters in the help message

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
    command = "nextflow run ${workflow.manifest.name} \
                    -profile conda \
                    --input </path/to/input_table> \
                    --output </path/to/output_folder> \
                    --reference </path/to/reference/folder>"

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
    // Create channel from input file provided through params.input
    //
    channel
    .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
    .map {
        meta, fastq_1, fastq_2 ->
            if (!fastq_2) {
                return [ meta.id, meta + [ paired_end:false ], [fastq_1]  ]
            } else {
                return [ meta.id, meta + [ paired_end:true ], [fastq_1, fastq_2]  ]
            }
    }
    .groupTuple()
    .map { samplesheet ->
        def (_id, meta, fastqs) = samplesheet
        workflow.profile.contains('test') ?
            [ meta[0], [ fastqs[0], fastqs[1] ].flatten() ] :
            validateInputSamplesheet(samplesheet)
    }
    .map { meta, fastqs ->
        [ meta, fastqs[0], fastqs[1] ]
    }
    .set { ch_samplesheet }


    //
    // Validate fusion tools provided by the user
    //
    fusiontools =  validateFusionTools(fusion_tools, model_pred)

    //
    // Build reference file channels
    //
    def reference_dir = reference.toString().replaceFirst(/\/$/, '')

    // configure test profile
    if (workflow.profile.contains('test')) {

        ref_fasta = "${reference_dir}/minigenome.fa"
        ref_gtf   = "${reference_dir}/minigenome.gtf"
        ref_tsl   = "${reference_dir}/minigenome.gtf.tsl"
        annot_db  = "${reference_dir}/minigenome.gff3.db"

        ch_ref_fasta = channel.value(file(ref_fasta, checkIfExists: true))
        ch_ref_gtf   = channel.value(file(ref_gtf, checkIfExists: true))
        ch_ref_tsl   = channel.value(file(ref_tsl, checkIfExists: true))
        ch_annot_db  = channel.value(file(annot_db, checkIfExists: true))

        stararriba_idx = file(
            "${reference_dir}/star_index.tar.gz",
            checkIfExists: true
        )

        starfusion_idx = file(
            "${reference_dir}/starfusion_index.tar.gz",
            checkIfExists: true
        )

        fusioncatcher_idx = file(
            "${reference_dir}/fusioncatcher_index.tar.gz",
            checkIfExists: true
        )

        UNTAR_STAR_INDEX(
            [[id:'test_idx'], stararriba_idx]
        )
        ch_stararriba_index = UNTAR_STAR_INDEX.out.untar
            .map { _meta, idx_path -> idx_path }

        UNTAR_STARFUSION_INDEX(
            [[id:'test_idx'], starfusion_idx]
        )
        ch_starfusion_index = UNTAR_STARFUSION_INDEX.out.untar
            .map { _meta, idx_path -> idx_path }

        UNTAR_FUSIONCATCHER_INDEX(
            [[ id:'test_idx'], fusioncatcher_idx]
        )
        ch_fusioncatcher_index = UNTAR_FUSIONCATCHER_INDEX.out.untar
            .map { _meta, idx_path -> idx_path }
    }
    else {

        ch_ref_fasta = channel.fromPath("${reference_dir}/*.dna.primary_assembly.fa", checkIfExists: true)
        ch_ref_gtf   = channel.fromPath("${reference_dir}/*.gtf", checkIfExists: true)
        ch_ref_tsl   = channel.fromPath("${reference_dir}/*.gtf.tsl", checkIfExists: true)
        ch_annot_db  = channel.fromPath("${reference_dir}/*.gff3.db", checkIfExists: true)

        ch_stararriba_index    = channel.value(file("${reference_dir}/star_index", checkIfExists: true))
        ch_starfusion_index    = channel.value(file("${reference_dir}/starfusion_index", checkIfExists: true))
        ch_fusioncatcher_index = channel.value(file("${reference_dir}/fusioncatcher_index", checkIfExists: true))
    }

    ch_model_threshold  = channel.value(model_threshold)
    ch_prediction_model = channel.value(
        file(
            "${projectDir}/assets/data/model/${model_prefix}.${model_pred}.rds",
            checkIfExists: true
        )
    )

    emit:

    samplesheet             = ch_samplesheet
    fusiontools             = fusiontools
    reference_fasta         = ch_ref_fasta
    reference_gtf           = ch_ref_gtf
    reference_tsl           = ch_ref_tsl
    annotation_db           = ch_annot_db
    starfusion_index        = ch_starfusion_index
    fusioncatcher_index     = ch_fusioncatcher_index
    stararriba_index        = ch_stararriba_index
    prediction_model        = ch_prediction_model
    model_threshold         = ch_model_threshold

    versions                = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Validate if the user provided fusion tools are a valid choice
//
def validateFusionTools(fusion_tools, model) {

    def valid_tools  = ['arriba', 'starfusion', 'fusioncatcher']
    def valid_models = ['EF_requant', 'EF_requant_type', 'EF_full']

    // Model corresponding to each validated fusion-tool combination.
    def tool_model_map = [
        'arriba,starfusion'               : ['EF_requant_type', 'EF_requant'],
        'arriba,fusioncatcher,starfusion' : 'EF_full',
    ]

    // Validate fusion tools
    if (!fusion_tools) {
        error "Please provide --fusion_tools. Valid options: ${valid_tools.join(', ')}"
    }

    def tools = fusion_tools
        .split(',')
        .collect { tool -> tool.trim().toLowerCase() }
        .findAll { tool -> tool }
        .unique()

    def invalid = tools - valid_tools

    if (invalid) {
        error "Invalid fusion tool(s): ${invalid.join(', ')}. " +
              "Valid options: ${valid_tools.join(', ')}"
    }

    // Make tool order irrelevant
    def tool_key = tools.unique().sort().join(',')

    // Automatically select model when --model is not provided
    if (!model) {

        if (!(tool_key in tool_model_map)) {
            error "No validated model is available for fusion tool combination: " +
                  "${tools.join(', ')}"
        }

        model = tool_model_map[tool_key]
    }

    // Validate model
    if (!(model in valid_models)) {
        error "Invalid model '${model}'. Valid options: ${valid_models.join(', ')}"
    }

    // Validate explicitly supplied model against tool combination
    def expected_model = tool_model_map[tool_key]

    if (expected_model && !(model in expected_model)) {
        error "Model '${model}' is not compatible with combination of fusion tools: " +
              "${tools.join(', ')}. Expected model: '${expected_model.join(', ')}'."
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

    def paired_end_values = metas.collect { meta -> meta.paired_end }.unique()

    if (paired_end_values.size() != 1) {
        error(
            "Please check input samplesheet -> " +
            "Multiple runs of a sample must be of the same datatype " +
            "i.e. single-end or paired-end: ${metas[0].id}"
        )
    }

    return [
        metas[0],
        fastqs.flatten()
    ]
}
