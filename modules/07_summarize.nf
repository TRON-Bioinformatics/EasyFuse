
process MERGE_DATA {
    tag "${name}"
    label 'process_single'
    publishDir "${params.output}/${name}", mode: 'copy'
    
    conda ("${baseDir}/environments/merging.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python:3.8.0--5e0e57f6a223cdda' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/796c2d0f22be8bd23bf91dd74661318db5545bc15c4ec63fd1207c95b2d5d22c/data' }"

    input:
      tuple val(name), path(detected_fusions), path(annot_fusions_csv), path(annot_fusions_csv_debug), path(annot_fusions_fasta), path(counts), path(read_stats)

    output:
      tuple val("${name}"), path("fusions.csv"), emit: merged_results

    script:
    """
    merge_data.py \
        --detected_fusions ${detected_fusions} \
        --context_seqs ${annot_fusions_csv} \
        --requant_counts ${counts} \
        --read_stats ${read_stats} \
        -o fusions.csv \
        --fusion_tools fusioncatcher,starfusion,arriba
    """
}

process PREDICTION {
    tag "${name}"
    label 'process_single'
    publishDir "${params.output}/${name}", mode: 'copy'

    conda ("${baseDir}/environments/prediction.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/r-base_r-optparse_r-randomforest_r-tidyverse:b6c2859feb55cb13' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fde5e222dd64b9eadebe8ceb397a8b2f3ed92afc55d7f5ff85824eceece7b7e9/data' }"

    input:
      tuple val(name), path(merged_results)

    output:
      tuple val("${name}"), path("fusions.pass.csv"), emit: predictions

    script:
    """
    R_model_prediction.R \
      --fusion_summary ${merged_results} \
      --model_file ${params.model_pred} \
      --prediction_threshold ${params.model_threshold} \
      --output fusions.pass.csv
    """
}