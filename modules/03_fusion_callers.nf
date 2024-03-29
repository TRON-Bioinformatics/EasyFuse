
process FUSION_CATCHER_INDEX {
    cpus 6
    memory "32g"
    tag "${name}"

    conda (params.enable_conda ? "${baseDir}/environments/fusioncatcher.yml" : null)

    script:
    """
    download-human-db.sh
    """
}


process FUSION_CATCHER {
    cpus 6
    memory "32g"
    tag "${name}"

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
    cpus 6
    memory "32g"
    tag "${name}"

    conda (params.enable_conda ? "${baseDir}/environments/starfusion.yml" : null)

    input:
      tuple val(name), path(chimeric_reads)

    output:
      tuple val("${name}"), path("star-fusion.fusion_predictions.tsv"), emit: fusions

    script:
    """
    STAR-Fusion \
        --chimeric_junction ${chimeric_reads} \
        --genome_lib_dir ${params.starfusion_index} \
        --CPU ${task.cpus} \
        --output_dir .
    """
}