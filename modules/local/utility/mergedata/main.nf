process MERGE_DATA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python:3.8.0--5e0e57f6a223cdda' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/796c2d0f22be8bd23bf91dd74661318db5545bc15c4ec63fd1207c95b2d5d22c/data' }"

    input:
    tuple val(meta),
          path(detected_fusions),
          path(annot_fusions_csv),
          path(annot_fusions_csv_debug),
          path(annot_fusions_fasta),
          path(counts),
          path(read_stats)

    output:
    tuple val(meta), path("fusions.csv"), emit: merged_results
    path("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    merge_data.py \\
        --detected_fusions ${detected_fusions} \\
        --context_seqs ${annot_fusions_csv} \\
        --requant_counts ${counts} \\
        --read_stats ${read_stats} \\
        -o fusions.csv \\
        --fusion_tools fusioncatcher,starfusion,arriba \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_data.py: \$(merge_data.py --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch fusions.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_data.py: \$(merge_data.py --version 2>&1)
    END_VERSIONS
    """
}
