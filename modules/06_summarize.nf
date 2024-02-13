
process MERGE_DATA {
    cpus 1
    memory "10g"
    tag "${name}"
    publishDir "${params.output}/${name}", mode: 'copy'
    
    conda (params.enable_conda ? "${baseDir}/environments/merging.yml" : null)

    input:
      tuple val(name), path(detected_fusions), path(annot_fusions_csv), path(annot_fusions_csv_debug), path(annot_fusions_fasta), path(counts), path(read_stats)

    output:
      tuple val("${name}"), path("fusions.csv"), emit: merged_results

    script:
    """
    merge_data.py \
        --input-fusions ${detected_fusions} \
        --input-fusion-context-seqs ${annot_fusions_csv} \
        --input-requant-counts ${counts} \
        --input-read-stats ${read_stats} \
        -o fusions.csv \
        --fusion-tools fusioncatcher,star,arriba
    """
}

process PREDICTION {
    cpus 1
    memory "10g"
    tag "${name}"
    publishDir "${params.output}/${name}", mode: 'copy'

    conda (params.enable_conda ? "${baseDir}/environments/prediction.yml" : null)

    input:
      tuple val(name), path(merged_results)

    output:
      tuple val("${name}"), path("fusions.pass.csv"), emit: predictions

    script:
    """
    R_model_prediction.R \
      --fusion_summary ${merged_results} \
      --model_file ${params.model_pred} \
      --prediction_threshold 0.5 \
      --output fusions.pass.csv
    """
}