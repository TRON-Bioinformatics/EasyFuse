process BPQUANT_CSV2FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bowtie2_bwa_pysam_samtools_pruned:dbf6a7df7fd19e94' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6cff3f1cbfc4d55dcafb0096f37bd334a48d175064b94e5d3d75245b29db43b/data'}"

    input:
    tuple val(meta), path(formatted_csv)

    output:
    tuple val(meta), path("${prefix}_seq_table.fasta"), emit: formatted_fasta
    path("versions.yml")                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bp_quant \\
        csv2fasta \\
        --input_csv ${formatted_csv} \\
	    --output_fasta ${prefix}_seq_table.fasta \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bp_quant: \$(bp_quant -h | grep 'version' | cut -d ' ' -f3)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_seq_table.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bp_quant: \$(bp_quant -h | grep 'version' | cut -d ' ' -f3)
    END_VERSIONS
    """
}
