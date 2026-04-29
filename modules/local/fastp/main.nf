process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda ("${moduleDir}/environment.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/fastp:0.23.4--f8cefc1e5f7a782e' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3fcff4f02e7e012e4bab124d64a2a50817dd64303998170127c8cf9c1968e10a/data' }"

    input:
    tuple val(meta), path(fastq1), path(fastq2)

    output:
    tuple val(meta),
    path("${prefix}_trimmed_R1.fastq.gz"),
    path("${prefix}_trimmed_R2.fastq.gz"), emit: trimmed_fastq
    path("versions.yml")                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def in1 = fastq1 ? "-in1 ${fastq1}" : ''
    def in2 = fastq2 ? "-in2 ${fastq2}" : ''
    def out1 = fastq1 ? "-out1 ${prefix}_trimmed_R1.fastq.gz" : ''
    def out2 = fastq2 ? "-out2 ${prefix}_trimmed_R2.fastq.gz" : ''

    """
    fastp \\
        ${in1} \\
        ${in2} \\
        ${out1} \\
        ${out2} \\
        --thread ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | cut -d ' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_trimmed_{R1,R2}.fastq
    gzip ${prefix}_trimmed_{R1,R2}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | cut -d ' ' -f2)
    END_VERSIONS
    """
}
