
process SUMMARIZE_DATA {
    cpus 2
    memory "10g"
    tag "${name}"
    publishDir "${params.output}/${name}", mode: 'copy'

    conda (params.enable_conda ? "${baseDir}/environments/easyfuse_src.yml" : null)

    input:
      tuple val(name), path(fusions), path(annot_csv), path(annot_fasta), path(cpms), path(counts), path(read_stats)

    output:
      tuple val("${name}"), path("fusions.csv"), path("fusions.pass.csv"), emit: predictions

    script:
    """
    easy-fuse summarize-data \
        --input-fusions ${fusions} \
        --input-fusion-context-seqs ${annot_csv} \
        --input-requant-cpm ${cpms} \
        --input-requant-counts ${counts} \
        --input-reads-stats ${read_stats} \
	--input-model-pred ${params.model_pred} \
        --model_predictions \
        -o . \
        --requant-mode best \
        --model-pred-threshold 0.5 \
        --fusion-tools fusioncatcher,star
    """
}
