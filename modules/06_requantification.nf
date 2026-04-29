process FUSION_FILTER {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/filtering.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pysam:0.22.0--a94c5bab35035aad' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64682c99fc92227f78f81a53c8d739b16e8b712c6c75a8909f159405cb29dbe1/data' }"

    input:
      tuple val(name), path(bam), path(annot_fusions_csv), path(annot_fusions_csv_debug), path(annot_fusions_fasta), path(read_stats)

    output:
      tuple val("${name}"), path("${name}.requantified.bam"), emit: bams

    script:
    """
    read_selection.py \
        --input ${bam} \
        --input2 ${annot_fusions_csv_debug} \
        --input_read_stats ${read_stats} \
        --output ${name}.requantified.bam
    """
}

process FUSION2CSV {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/filtering.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pysam:0.22.0--a94c5bab35035aad' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64682c99fc92227f78f81a53c8d739b16e8b712c6c75a8909f159405cb29dbe1/data' }"

    input:
      tuple val(name), path(annot_fusions_csv), path(annot_fusions_csv_debug), path(annot_fusions_fasta)

    output:
      tuple val("${name}"), path("seq_table.csv"), emit: formatted_csv

    script:
    """
    format_seq_table.py \
        --input_table ${annot_fusions_csv_debug} \
        --output_table seq_table.csv
    """
}

process CSV2FASTA {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/requantification.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bowtie2_bwa_pysam_samtools_pruned:dbf6a7df7fd19e94' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6cff3f1cbfc4d55dcafb0096f37bd334a48d175064b94e5d3d75245b29db43b/data'}"

    input:
      tuple val(name), path(formatted_csv)

    output:
      tuple val("${name}"), path("seq_table.fasta"), emit: formatted_fasta

    script:
    """
    bp_quant csv2fasta \
        --input_csv ${formatted_csv} \
	      --output_fasta seq_table.fasta
    """

}


process STAR_INDEX {
    tag "${name}"
    label 'process_low'

    conda ("${baseDir}/environments/requantification.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bowtie2_bwa_pysam_samtools_pruned:dbf6a7df7fd19e94' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6cff3f1cbfc4d55dcafb0096f37bd334a48d175064b94e5d3d75245b29db43b/data'}"
  
    input:
      tuple val(name), path(formatted_fasta)

    output:
      tuple val("${name}"), path("star_index"), emit: star_index

    script:
    """
    mkdir star_index
    bp_quant index \
        -i ${formatted_fasta} \
	      -o star_index \
	      -t 2 \
	      -m star
    """

}

process STAR_CUSTOM {
    tag "${name}"
    label 'process_medium'

    conda ("${baseDir}/environments/requantification.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bowtie2_bwa_pysam_samtools_pruned:dbf6a7df7fd19e94' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6cff3f1cbfc4d55dcafb0096f37bd334a48d175064b94e5d3d75245b29db43b/data'}"

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
      -t ${task.cpus} \
      -m star

    mv Aligned.out.sam ${name}.sam
    """
}

process READ_COUNT {
  tag "${name}"
  label 'process_single'

  conda ("${baseDir}/environments/requantification.yml")
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bowtie2_bwa_pysam_samtools_pruned:dbf6a7df7fd19e94' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6cff3f1cbfc4d55dcafb0096f37bd334a48d175064b94e5d3d75245b29db43b/data'}"

  input:
    tuple val(name), path(bam), path(formatted_csv)

  output:
    tuple val("${name}"), path("quantification.tsv"), emit: counts

  script:
  """
  mkdir results
  bp_quant count \
      --input ${bam} \
      --seq_table ${formatted_csv} \
      --output . \
      --bp_distance 10 \
      --allow_mismatches
  """
}