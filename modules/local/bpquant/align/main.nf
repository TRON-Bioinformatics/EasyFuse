process BPQUANT_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bowtie2_bwa_pysam_samtools_pruned:dbf6a7df7fd19e94' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6cff3f1cbfc4d55dcafb0096f37bd334a48d175064b94e5d3d75245b29db43b/data'}"

    input:
    tuple val(meta), path(fastq1), path(fastq2), path(star_index, stageAs: "star_index/")

    output:
    tuple val(meta), path("${prefix}.sam")          , emit: bam
    tuple val(meta), path("${prefix}.Log.final.out"), emit: read_stats
    path("versions.yml")                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bp_quant \\
        align \\
        -1 ${fastq1} \\
        -2 ${fastq2} \\
        -i ${star_index} \\
        -o . \\
        -t ${task.cpus} \\
        -m star \\
        ${args}

    mv Aligned.out.sam ${prefix}.sam
    mv Log.final.out ${prefix}.Log.final.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bpquant: \$(bp_quant -h | grep 'version' | cut -d ' ' -f3)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.sam
    touch ${prefix}.Log.final.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bpquant: \$(bp_quant -h | grep 'version' | cut -d ' ' -f3)
    END_VERSIONS
    """
}
