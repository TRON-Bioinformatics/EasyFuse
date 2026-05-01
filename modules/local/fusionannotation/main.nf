process FUSIONANNOTATER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/gffutils_biopython_python-xxhash_python:458dcfe93321068c' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/30840bab8184dfe45a8e2d0f41c869e5fcc447c4e83f0c9db5e74e9e3c07fa87/data' }"

    input:
    tuple val(meta), path(fusions), path(annotation_db)

    output:
    tuple val(meta),
          path("${prefix}_annotated_fusions.csv"),
          path("*annotated_fusions.csv.debug"),
          path("*annotated_fusions.csv.fasta"), emit: annot_fusions
    path("versions.yml")                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    fusionannotator.py \\
        --detected_fusions ${fusions} \\
        --annotation_db ${annotation_db} \\
        --out_csv ${prefix}_annotated_fusions.csv \\
        --genome_fasta ${params.fasta} \\
        --tsl_info ${params.reference_tsl} \\
        --cis_near_dist 1000000 \\
        --context_seq_len 400 \\
        --tsl_filter_level 4,5,NA \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusionannotation: \$(fusionannotator.py --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_annotated_fusions.csv
    touch ${prefix}_annotated_fusions.csv.debug
    touch ${prefix}_annotated_fusions.csv.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusionannotation: \$(fusionannotator.py --version)
    END_VERSIONS
    """
}
