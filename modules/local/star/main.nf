process STAR_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda ("${moduleDir}/environment.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9b/9b8ecb2f9a77b5e7573ef6fae2f4c2e771064f7a129ed1329913c1025c33f365/data' }"

    input:
    tuple val(meta), path(fastq1), path(fastq2), path(star_index, stageAs: "star_index/")

    output:
    tuple val(meta), path("${prefix}.bam")                  , emit: bam
    tuple val(meta), path("${prefix}.Log.final.out")        , emit: read_stats
    tuple val(meta), path("${prefix}.Chimeric.out.junction"), emit: chimeric_reads
    path("versions.yml")                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    STAR \\
        --genomeDir ${star_index} \\
        --outFileNamePrefix ${prefix}. \\
        --readFilesCommand zcat \\
        --readFilesIn ${fastq1} ${fastq2} \\
        --outFilterMultimapNmax 1000 \\
        --outSAMmultNmax 1 \\
        --chimSegmentMin 10 \\
        --chimJunctionOverhangMin 10 \\
        --chimOutJunctionFormat 1 \\
        --alignSJDBoverhangMin 10 \\
        --alignMatesGapMax 200000 \\
        --alignIntronMax 200000 \\
        --chimSegmentReadGapMax 3 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --seedSearchStartLmax 20 \\
        --winAnchorMultimapNmax 50 \\
        --outSAMtype BAM Unsorted \\
        --chimOutType Junctions WithinBAM \\
        --outSAMunmapped Within KeepPairs \\
        --runThreadN ${task.cpus} \\
        ${args}

    mv ${prefix}.Aligned.out.bam ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(star --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Chimeric.out.junction

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(star --version)
    END_VERSIONS
    """
}
