process FUSIONCATCHER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fusioncatcher:1.33--hdfd78af_4':
        'biocontainers/fusioncatcher:1.33--hdfd78af_4' }"

    input:
    tuple val(meta), path(fastq1), path(fastq2), path(fusioncatcher_index, stageAs: "fusioncatcher_index/")

    output:
    tuple val(meta),
    path("${prefix}/summary_candidate_fusions.txt"),
    path("${prefix}/final-list_candidate-fusion-genes.txt"), emit: fusions
    path("versions.yml")                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def run_type    = fastq2 ? "" : "--single-end"
    def input_fastq = fastq2 ? "${fastq1},${fastq2}" : "${fastq1}"

    """
    fusioncatcher \\
        --data ${fusioncatcher_index} \\
        --input ${input_fastq} \\
        --output ${prefix} \\
        ${run_type} \\
        -p ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(fusioncatcher --version | cut -d ' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/summary_candidate_fusions.txt
    touch ${prefix}/final-list_candidate-fusion-genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(fusioncatcher --version | cut -d ' ' -f2)
    END_VERSIONS
    """
}
