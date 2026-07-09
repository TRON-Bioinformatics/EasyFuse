process FUSION_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pysam:0.22.0--a94c5bab35035aad' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64682c99fc92227f78f81a53c8d739b16e8b712c6c75a8909f159405cb29dbe1/data' }"

    input:
    tuple val(meta), path(bam), path(annot_fusions_csv), path(annot_fusions_csv_debug), path(annot_fusions_fasta), path(read_stats)

    output:
    tuple val(meta), path("${prefix}.requantified.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    read_selection.py \\
        --input ${bam} \\
        --input2 ${annot_fusions_csv_debug} \\
        --input_read_stats ${read_stats} \\
        --output ${prefix}.requantified.bam \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.requantified.bam
    """
}
