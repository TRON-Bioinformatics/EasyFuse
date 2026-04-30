process ARRIBA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/arriba:2.4.0--9680480f3735ac7f' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fb/fbbd3ccedb1663939f2ca075a071e75b0d1c60f19a4cd46dd9ffe371f133105a/data' }"

    input:
    tuple val(meta), path(bam), path(gtf), path(fasta)

    output:
    tuple val(meta), path("${prefix}_fusions.tsv"), emit: fusions
    path("versions.yml")                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    arriba \\
        -x ${bam} \\
        -g ${gtf} \\
        -a ${fasta} \\
        -o ${prefix}_fusions.tsv \\
        -O ${prefix}_fusions.discarded.tsv \\
        -f blacklist \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba | head -1 | cut -d ' ' -f4)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_fusions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba | head -1 | cut -d ' ' -f4)
    END_VERSIONS
    """
}
