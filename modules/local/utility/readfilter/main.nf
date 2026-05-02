process READ_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pysam:0.22.0--a94c5bab35035aad' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64682c99fc92227f78f81a53c8d739b16e8b712c6c75a8909f159405cb29dbe1/data' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.filtered.bam"), emit: filtered_bam
    path("versions.yml")                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    fusionreadfilter.py \\
        --input ${bam} \\
        --output ${prefix}.filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusionreadfilter.py: \$(fusionreadfilter.py --version 2>&1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusionreadfilter.py: \$(fusionreadfilter.py --version 2>&1)
    END_VERSIONS
    """
}
