process FUSION2CSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pysam:0.22.0--a94c5bab35035aad' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64682c99fc92227f78f81a53c8d739b16e8b712c6c75a8909f159405cb29dbe1/data' }"

    input:
    tuple val(meta), path(annot_fusions_csv), path(annot_fusions_csv_debug), path(annot_fusions_fasta)

    output:
    tuple val(meta), path("${prefix}_seq_table.csv"), emit: formatted_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    format_seq_table.py \\
        --input_table ${annot_fusions_csv_debug} \\
        --output_table ${prefix}_seq_table.csv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_seq_table.csv
    """
}
