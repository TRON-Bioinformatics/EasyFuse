
process STAR {
    tag "${name}"
    label 'process_medium'

    conda ("${baseDir}/environments/alignment.yml")

    input:
      tuple val(name), path(fastq1), file(fastq2)

    output:
      tuple val("${name}"), path("${name}.bam"), emit: bams
      tuple val("${name}"), path("${name}.Chimeric.out.junction"), emit: chimeric_reads
      tuple val("${name}"), path("${name}.Log.final.out"), emit: read_stats

    script:
    """
    STAR --genomeDir ${params.star_index} \
        --outFileNamePrefix ${name}. \
        --readFilesCommand zcat \
        --readFilesIn ${fastq1} ${fastq2} \
        --outFilterMultimapNmax 1000 \
        --outSAMmultNmax 1 \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --chimOutJunctionFormat 1 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 200000 \
        --alignIntronMax 200000 \
        --chimSegmentReadGapMax 3 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --seedSearchStartLmax 20 \
        --winAnchorMultimapNmax 50 \
        --outSAMtype BAM Unsorted \
        --chimOutType Junctions WithinBAM \
        --outSAMunmapped Within KeepPairs \
        --runThreadN ${task.cpus}

    mv ${name}.Aligned.out.bam ${name}.bam
    """
}


process STAR_ARRIBA {
    tag "${name}"
    label 'process_medium'

    conda ("${baseDir}/environments/alignment.yml")

    input:
      tuple val(name), path(fastq1), file(fastq2)

    output:
      tuple val("${name}"), path("${name}.bam"), emit: bams

    script:
    """
    STAR --genomeDir ${params.star_index} \
        --outFileNamePrefix ${name}. \
        --readFilesCommand zcat \
        --readFilesIn ${fastq1} ${fastq2} \
        --outFilterMultimapNmax 50 \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 50 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMtype BAM Unsorted \
        --chimOutType WithinBAM HardClip \
        --outSAMunmapped Within \
        --runThreadN ${task.cpus}

    mv ${name}.Aligned.out.bam ${name}.bam
    """
}


process READ_FILTER {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/filtering.yml")

    input:
      tuple val(name), path(bam)

    output:
      tuple val("${name}"), path("${name}.filtered.bam"), emit: bams

    script:
    """
    fusionreadfilter.py \
    --input ${bam} \
    --output ${name}.filtered.bam
    """
}

process BAM2FASTQ {
    tag "${name}"
    label 'process_medium'

    conda ("${baseDir}/environments/samtools.yml")

    input:
      tuple val(name), path(bam)

    output:
      tuple val("${name}"), path("${name}.read1.fastq.gz"), path("${name}.read2.fastq.gz"), emit: fastqs

    script:
    """
    samtools fastq \
    -0 ${name}.other.fastq.gz \
    -1 ${name}.read1.fastq.gz \
    -2 ${name}.read2.fastq.gz \
    --threads ${task.cpus} \
    ${bam}
    """
}