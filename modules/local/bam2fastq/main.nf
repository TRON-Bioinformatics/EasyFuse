process BAM2FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/samtools:1.9--ff876e25d460de68' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/597f0179350adad0f4fa6b045f64f469fc04bf6c1a26c870bf900f2832bf260f/data' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.read1.fastq.gz"), path("${prefix}.read2.fastq.gz"), emit: fastqs
    path("versions.yml")                                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools \\
        fastq \\
        -0 ${prefix}.other.fastq.gz \\
        -1 ${prefix}.read1.fastq.gz \\
        -2 ${prefix}.read2.fastq.gz \\
        --threads ${task.cpus} \\
        ${bam} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2fastq (samtools): \$(samtools version | sed '1!d;s/.* //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.read1.fastq ${prefix}.read2.fastq
    gzip ${prefix}.read1.fastq ${prefix}.read2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2fastq (samtools): \$(samtools version | sed '1!d;s/.* //')
    END_VERSIONS
    """
}
