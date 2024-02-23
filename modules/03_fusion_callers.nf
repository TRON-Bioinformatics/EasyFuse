
process FUSION_CATCHER_INDEX {
    tag "${name}"
    label 'process_low'

    conda (params.enable_conda ? "${baseDir}/environments/fusioncatcher.yml" : null)

    script:
    """
    download-human-db.sh
    """
}


process FUSION_CATCHER {
    tag "${name}"
    label 'process_medium'

    conda (params.enable_conda ? "${baseDir}/environments/fusioncatcher.yml" : null)

    input:
      tuple val(name), path(fastq1), file(fastq2)

    output:
      tuple val("${name}"), path("summary_candidate_fusions.txt"), path("final-list_candidate-fusion-genes.txt"), emit: fusions

    script:
    """
    
    fusioncatcher \
        --data ${params.fusioncatcher_index} \
        --input ${fastq1},${fastq2} \
        --output . \
        -p ${task.cpus}
    """
}

process STAR_FUSION {
    tag "${name}"
    label 'process_medium'

    conda (params.enable_conda ? "${baseDir}/environments/starfusion.yml" : null)

    input:
      tuple val(name), path(fastq1), file(fastq2)

    output:
      tuple val("${name}"), path("star-fusion.fusion_predictions.tsv"), emit: fusions

    script:
    """
    STAR-Fusion \
        --left_fq ${fastq1} \
        --right_fq ${fastq2} \
        --genome_lib_dir ${params.starfusion_index} \
        --CPU ${task.cpus} \
        --output_dir .
    """
}

process ARRIBA {
    tag "${name}"
    label 'process_low'

    conda (params.enable_conda ? "${baseDir}/environments/arriba.yml" : null)

    input:
      tuple val(name), path(bam)

    output:
      tuple val("${name}"), path("fusions.tsv"), path("fusions.discarded.tsv"), emit: fusions

    script:
    """
    arriba \
        -x ${bam} \
        -g ${params.gtf} \
        -a ${params.fasta} \
        -o fusions.tsv \
        -O fusions.discarded.tsv \
        -f blacklist
    """
}