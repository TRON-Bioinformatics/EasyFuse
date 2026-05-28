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
    val(fusion_tools)

    output:
    tuple val(meta), path("fusions.csv"), emit: merged_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tool_list = []
    if (fusion_tools.run_fusioncatcher) tool_list.add('fusioncatcher')
    if (fusion_tools.run_starfusion)    tool_list.add('starfusion')
    if (fusion_tools.run_arriba)        tool_list.add('arriba')

    """
    merge_data.py \\
        --detected_fusions ${detected_fusions} \\
        --context_seqs ${annot_fusions_csv} \\
        --requant_counts ${counts} \\
        --read_stats ${read_stats} \\
        -o fusions.csv \\
        --fusion_tools ${tool_list.join(',')} \\
        ${args}
    """

    stub:

    """
    touch fusions.csv
    """
}
