process PARSE_FUSIONCATCHER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python:3.8.0--5e0e57f6a223cdda' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/796c2d0f22be8bd23bf91dd74661318db5545bc15c4ec63fd1207c95b2d5d22c/data' }"

    input:
    tuple val(meta), path(fusion_catcher_1), path(fusion_catcher_2)

    output:
    tuple val(meta), path("${prefix}_fusioncatcher.csv"), emit: fusions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    parse_tool.py \\
        --input_file ${fusion_catcher_1} \\
        --input_file2 ${fusion_catcher_2} \\
        --output_file ${prefix}_fusioncatcher.csv \\
        --tool fusioncatcher \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_fusioncatcher.csv
    """
}
