process FUSION_PARSER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python:3.8.0--5e0e57f6a223cdda' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/796c2d0f22be8bd23bf91dd74661318db5545bc15c4ec63fd1207c95b2d5d22c/data' }"

    input:
    tuple val(meta), path(fusioncatcher_fusions), path(starfusion_fusions), path(arriba_fusions)

    output:
    tuple val(meta), path("${prefix}_Detected_Fusions.csv"), emit: fusions
    path("versions.yml")                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def arriba        = arriba_fusions        ? "--tool arriba ${arriba_fusions}" : ''
    def starfusion    = starfusion_fusions    ? "--tool starfusion ${starfusion_fusions}" : ''
    def fusioncatcher = fusioncatcher_fusions ? "--tool fusioncatcher ${fusioncatcher_fusions}" : ''

    """
    fusiontoolparser.py \\
        ${fusioncatcher} \\
	    ${starfusion} \\
        ${arriba} \\
        --output ${prefix}_Detected_Fusions.csv \\
        --sample ${prefix} \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsefusions: \$(fusiontoolparser.py --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_Detected_Fusions.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsefusions: \$(fusiontoolparser.py --version)
    END_VERSIONS
    """
}
