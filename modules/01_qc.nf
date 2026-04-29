process FASTP {
    tag "${name}"
    label 'process_low'

    conda ("${baseDir}/environments/qc.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/fastp:0.23.4--f8cefc1e5f7a782e' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3fcff4f02e7e012e4bab124d64a2a50817dd64303998170127c8cf9c1968e10a/data' }"

    input:
      tuple val(name), path(fastq1), path(fastq2)

    output:
      tuple val("${name}"), path("trimmed_R1.fastq.gz"), path("trimmed_R2.fastq.gz"), emit: trimmed_fastq

    script:
    """
        fastp \\
            -i ${fastq1} \\
            -I ${fastq2} \\
            -o trimmed_R1.fastq.gz \\
            -O trimmed_R2.fastq.gz \\
            --thread ${task.cpus}
    """
}
