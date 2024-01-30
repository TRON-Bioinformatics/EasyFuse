process FUSION_FILTER {
    cpus 1
    memory "8g"
    tag "${name}"

    conda (params.enable_conda ? "${baseDir}/environments/filtering.yml" : null)

    input:
      tuple val(name), path(bam), path(annot_fusions_1), path(annot_fusions_2), path(read_stats)

    output:
      tuple val("${name}"), path("${name}.requantified.bam"), emit: bams

    script:
    """
    read_selection.py \
        --input ${bam} \
        --input2 ${annot_fusions_1} \
        --input-reads-stats ${read_stats} \
        --output ${name}.requantified.bam
    """
}

process FUSION2CSV {
    cpus 1
    memory "1g"
    tag "${name}"

    input:
      tuple val(name), path(annot_fusions_1), path(annot_fusions_2)

    output:
      tuple val("${name}"), path("seq_table.csv"), emit: annot_csv

    script:
    """
    format_seq_table.py \
        --input_table ${annot_fusions_1} \
        --output_table seq_table.csv
    """
}

process CSV2FASTA {
    cpus 1
    memory "1g"
    tag "${name}"

    conda (params.enable_conda ? "${baseDir}/environments/requantification.yml" : null)

    input:
      tuple val(name), path(annot_csv)

    output:
      tuple val("${name}"), path("seq_table.fasta"), emit: annot_fasta

    script:
    """
    bp_quant csv2fasta \
        --input_csv ${annot_csv} \
	      --output_fasta seq_table.fasta
    """

}


process STAR_INDEX {
    cpus 2
    memory "8g"
    tag "${name}"

    conda (params.enable_conda ? "${baseDir}/environments/requantification.yml" : null)
    
    input:
      tuple val(name), path(annot_fusions_1), path(annot_fusions_2)

    output:
      tuple val("${name}"), path("star_index"), emit: star_index

    script:
    """
    mkdir star_index
    bp_quant index \
        -i ${annot_fusions_2} \
	      -o star_index \
	      -t 2 \
	      -m star
    """

}

process STAR_CUSTOM {
    cpus 6
    memory "32g"
    tag "${name}"

    conda (params.enable_conda ? "${baseDir}/environments/requantification.yml" : null)

    input:
      tuple val(name), path(fastq1), file(fastq2), path(star_index)

    output:
      tuple val("${name}"), path("${name}.sam"), emit: bams
      tuple val("${name}"), path("Log.final.out"), emit: read_stats

    script:
    """
    bp_quant align \
      -1 ${fastq1} \
      -2 ${fastq2} \
      -i ${star_index} \
      -o . \
      -t 6 \
      -m star

    mv Aligned.out.sam ${name}.sam
    """
}

process READ_COUNT {
  cpus 6
  memory "50g"
  tag "${name}"

  conda (params.enable_conda ? "${baseDir}/environments/requantification.yml" : null)

  input:
    tuple val(name), path(bam), path(annot_csv)

  output:
    tuple val("${name}"), path("quantification.tsv"), emit: counts

  script:
  """
  mkdir results
  bp_quant count \
      --input ${bam} \
      --seq_table ${annot_csv} \
      --output . \
      --bp_distance 10 \
      --allow_mismatches

  """
}