
process SUMMARIZE_DATA {
    cpus 2
    memory "10g"
    tag "${name}"
    publishDir "${params.output}/${name}", mode: 'copy'

    conda (params.enable_conda ? "environments/summarize_data.yml" : null)

    input:
      tuple val(name), path(fusions), path(annot_csv), path(annot_fasta), path(cpms), path(counts), path(read_stats)

    output:
      tuple val("${name}"), path("fusions_1.csv"), path("fusions.pass.csv"), path("fusions.pass.all.csv"), emit: predictions

    script:
    """
    easy-fuse summarize-data \
        --input-fusions ${fusions} \
        --input-fusion-context-seqs ${annot_csv} \
        --input-requant-cpm ${cpms} \
        --input-requant-counts ${counts} \
        --input-reads-stats ${read_stats} \
        --model_predictions \
        -o . \
        --requant-mode best \
        --model-pred-threshold 0.5 \
        --fusion-tools fusioncatcher,star
    """
}
