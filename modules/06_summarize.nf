
process SUMMARIZE_DATA {
    cpus 1
    memory "10g"
    tag "${name}"
    publishDir "${params.output}/${name}", mode: 'copy'

    conda (params.enable_conda ? "${baseDir}/environments/summarize.yml" : null)

    input:
      tuple val(name), path(fusions), path(annot_csv), path(annot_fasta), path(counts), path(read_stats)

    output:
      tuple val("${name}"), path("fusions.csv"), path("fusions.pass.csv"), emit: predictions

    script:
    """
    summarize_data.py \
        --input-fusions ${fusions} \
        --input-fusion-context-seqs ${annot_csv} \
        --input-requant-counts ${counts} \
        --input-read-stats ${read_stats} \
	      --input-model-pred ${params.model_pred} \
        --model_predictions \
        -o . \
        --model-pred-threshold 0.5 \
        --fusion-tools fusioncatcher,star,arriba
    """
}
