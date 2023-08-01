
process STAR_INDEX {
    cpus 2
    memory "8g"
    tag "${name}"

    // NOTE: we really need STAR and samtools here, these are not added as dependencies to the bioconda package
    conda (params.enable_conda ? "${baseDir}/environments/requantification.yml" : null)

    input:
      tuple val(name), path(annot_csv), path(annot_fasta)

    output:
      tuple val("${name}"), path("star_index"), emit: star_index

    script:

    """
    mkdir star_index
    easy-fuse star-index \
        -i ${annot_fasta} \
        -o star_index \
        -t 1 \
        -b STAR
    """
}

process FUSION_FILTER {
    cpus 2
    memory "8g"
    tag "${name}"

    conda (params.enable_conda ? "${baseDir}/environments/easyfuse_src.yml" : null)

    input:
      tuple val(name), path(bam), path(annot_fusions_1), path(annot_fusions_2), path(read_stats)

    output:
      tuple val("${name}"), path("${name}.requantified.bam"), emit: bams

    script:

    """
    easy-fuse requantify-filter \
        --input ${bam} \
        --input2 ${annot_fusions_1} \
        --input-reads-stats ${read_stats} \
        --output ${name}.requantified.bam
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
      tuple val("${name}"), path("${name}.bam"), path("${name}.bam.bai"), emit: bams
      tuple val("${name}"), path("${name}.Log.final.out"), emit: read_stats

    script:
    """

    STAR --genomeDir ${star_index} \
        --outFileNamePrefix ${name}. \
        --readFilesCommand zcat \
        --readFilesIn ${fastq1} ${fastq2} \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax -1 \
        --outSAMattributes Standard \
        --outSAMunmapped None \
        --outFilterMismatchNoverLmax 0.02 \
        --runThreadN ${task.cpus}

    mv ${name}.Aligned.sortedByCoord.out.bam ${name}.bam

    samtools index ${name}.bam
    """
}

process READ_COUNT {
  cpus 6
  memory "50g"
  tag "${name}"

  conda (params.enable_conda ? "${baseDir}/environments/easyfuse_src.yml" : null)

  input:
    tuple val(name), path(bam), path(bam_index), path(read_stats)

  output:
    tuple val("${name}"), path("requant.csv"), path("requant.csv.counts"), emit: counts

  script:
  """

  easy-fuse requantify \
      --input ${bam} \
      --output requant.csv \
      --bp_distance 10 \
      --input-reads-stats ${read_stats}

  """
}