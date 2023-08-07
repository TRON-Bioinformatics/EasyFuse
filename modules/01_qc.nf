
process FASTQC {
    cpus 2
    memory "8g"
    tag "${name}"

    conda (params.enable_conda ? "${baseDir}/environments/qc.yml" : null)

    input:
    tuple val(name),
        path(fastq1),
        path(fastq2)

    output:
    tuple val("${name}"),
        path("${fastq1.simpleName}_fastqc_data.txt"),
        path("${fastq2.simpleName}_fastqc_data.txt"), emit: qc_data

    """
    fastqc --nogroup --extract -t 6 -o . $fastq1 $fastq2
    mv ${fastq1.simpleName}_fastqc/fastqc_data.txt ${fastq1.simpleName}_fastqc_data.txt
    mv ${fastq2.simpleName}_fastqc/fastqc_data.txt ${fastq2.simpleName}_fastqc_data.txt
    """
}

process FASTQC_PARSER {
    cpus 2
    memory "8g"
    tag "${name}"

    conda (params.enable_conda && ! params.disable_pyeasyfuse_conda ? "${baseDir}/environments/easyfuse_src.yml" : null)

    input:
    tuple val(name),
        path(fastq1_qc_data),
        path(fastq2_qc_data)

    output:
    tuple val("${name}"),
        path("qc_table.txt"), emit: qc_table

    """
    easy-fuse qc-parser -i $fastq1_qc_data $fastq2_qc_data -o qc_table.txt
    """
}

process SKEWER {
    cpus 2
    memory "8g"
    tag "${name}"

    conda (params.enable_conda ? "${baseDir}/environments/easyfuse_src.yml" : null)

    input:
    tuple val(name),
        path(fastq1),
        path(fastq2),
        path(qc_table)

    output:
    tuple val("${name}"),
        path("out_file-trimmed-pair1.fastq.gz"),
        path("out_file-trimmed-pair2.fastq.gz"), emit: trimmed_fastq

    """
    easy-fuse skewer-wrapper -q ${qc_table} \
        -i ${fastq1} ${fastq2} \
        -o . \
        -b skewer \
        -m 0.75
    """
}
