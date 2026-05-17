process STARFUSION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/star-fusion:1.12.0--359bb9f50e24aa17' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a8566c3e2ecd2afdb44563a65e894975fb746f335df651e76bc6c127cc62029c/data' }"

    input:
    tuple val(meta), path(fastq1), path(fastq2), path(starfusion_index, stageAs: "starfusion_index/")

    output:
    tuple val(meta), path("${prefix}/star-fusion.fusion_predictions.tsv"), emit: fusions
    path("versions.yml")                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    STAR-Fusion \\
        --left_fq ${fastq1} \\
        --right_fq ${fastq2} \\
        --genome_lib_dir ${starfusion_index} \\
        --CPU ${task.cpus} \\
        --output_dir ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        starfusion: \$(STAR-Fusion --version | tr -d '\\n' | cut -d ' ' -f3)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/star-fusion.fusion_predictions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        starfusion: \$(STAR-Fusion --version | tr -d '\\n' | cut -d ' ' -f3)
    END_VERSIONS
    """
}
