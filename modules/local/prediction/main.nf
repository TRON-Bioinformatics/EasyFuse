process PREDICTION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/r-base_r-optparse_r-randomforest_r-tidyverse:b6c2859feb55cb13' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fde5e222dd64b9eadebe8ceb397a8b2f3ed92afc55d7f5ff85824eceece7b7e9/data' }"

    input:
    tuple val(meta), path(merged_results)

    output:
    tuple val(meta), path("fusions.pass.csv"), emit: predictions
    path("versions.yml")                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def _prefix = task.ext.prefix ?: "${meta.id}"

    """
    R_model_prediction.R \\
      --fusion_summary ${merged_results} \\
      --model_file ${params.model_pred} \\
      --prediction_threshold ${params.model_threshold} \\
      --output fusions.pass.csv \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R_model_prediction.R: \$(R_model_prediction.R --version)
    END_VERSIONS
    """

    stub:

    """
    touch fusions.pass.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R_model_prediction.R: \$(R_model_prediction.R --version)
    END_VERSIONS
    """
}
