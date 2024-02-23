
process FASTP {
    tag "${name}"
    label 'process_medium'

    conda (params.enable_conda ? "${baseDir}/environments/qc.yml" : null)

    input:
      tuple val(name), path(fastq1), path(fastq2)

    output:
      tuple val("${name}"), path("trimmed_R1.fastq.gz"), path("trimmed_R2.fastq.gz"), emit: trimmed_fastq

    """
    fastp \
        -i ${fastq1} \
        -I ${fastq2} \
        -o trimmed_R1.fastq.gz \
        -O trimmed_R2.fastq.gz \
        --thread 6
    """
}
