
process FUSION_CATCHER {
    tag "${name}"
    label 'process_medium'

    conda ("${baseDir}/environments/fusioncatcher.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/fusioncatcher:1.33--4733482b637ef92f':
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d53f36e9e01d14a0ae8e15f8046f52b2883c970c27fe43fdfbd9440a55f5403f/data' }"

    input:
      tuple val(name), path(fastq1), file(fastq2)
      path(fusioncatcher_index, stageAs: "fusioncatcher_index/")

    output:
      tuple val("${name}"), path("summary_candidate_fusions.txt"), path("final-list_candidate-fusion-genes.txt"), emit: fusions

    script:
    """
    
    fusioncatcher \
        --data ${fusioncatcher_index} \
        --input ${fastq1},${fastq2} \
        --output . \
        -p ${task.cpus}
    """
}

process STAR_FUSION {
    tag "${name}"
    label 'process_medium'

    conda ("${baseDir}/environments/starfusion.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/star-fusion:1.12.0--359bb9f50e24aa17' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a8566c3e2ecd2afdb44563a65e894975fb746f335df651e76bc6c127cc62029c/data' }"

    input:
      tuple val(name), path(fastq1), file(fastq2)
      path(starfusion_index, stageAs: "starfusion_index/")

    output:
      tuple val("${name}"), path("star-fusion.fusion_predictions.tsv"), emit: fusions

    script:
    """
    STAR-Fusion \
        --left_fq ${fastq1} \
        --right_fq ${fastq2} \
        --genome_lib_dir ${starfusion_index} \
        --CPU ${task.cpus} \
        --output_dir .
    """
}

process ARRIBA {
    tag "${name}"
    label 'process_medium'

    conda ("${baseDir}/environments/arriba.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/arriba:2.4.0--9680480f3735ac7f' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fb/fbbd3ccedb1663939f2ca075a071e75b0d1c60f19a4cd46dd9ffe371f133105a/data' }"

    input:
      tuple val(name), path(bam)
      path(gtf)
      path(fasta)

    output:
      tuple val("${name}"), path("fusions.tsv"), emit: fusions

    script:
    """
    arriba \
        -x ${bam} \
        -g ${gtf} \
        -a ${fasta} \
        -o fusions.tsv \
        -O fusions.discarded.tsv \
        -f blacklist
    """
}
