process BPQUANT_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bowtie2_bwa_pysam_samtools_pruned:dbf6a7df7fd19e94' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6cff3f1cbfc4d55dcafb0096f37bd334a48d175064b94e5d3d75245b29db43b/data'}"

    input:
    tuple val(meta), path(formatted_fasta)

    output:
    tuple val(meta), path("star_index"), emit: star_index
    path("versions.yml")               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p star_index

    bp_quant \\
        index \\
        -i ${formatted_fasta} \\
	    -o star_index \\
	    -t ${task.cpus} \\
	    -m star \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bpquant: \$(bp_quant -h | grep 'version' | cut -d ' ' -f3)
    END_VERSIONS
    """

    stub:

    """
    mkdir -p star_index/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bpquant: \$(bp_quant -h | grep 'version' | cut -d ' ' -f3)
    END_VERSIONS
    """
}
