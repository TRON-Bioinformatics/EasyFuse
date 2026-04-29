process STAR {
    tag "${name}"
    label 'process_medium'

    conda ("${baseDir}/environments/alignment.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9b/9b8ecb2f9a77b5e7573ef6fae2f4c2e771064f7a129ed1329913c1025c33f365/data' }"

    input:
      tuple val(name), path(fastq1), path(fastq2)
      path(star_index, stageAs: "star_index/")

    output:
      tuple val("${name}"), path("${name}.bam"), emit: bams
      tuple val("${name}"), path("${name}.Chimeric.out.junction"), emit: chimeric_reads
      tuple val("${name}"), path("${name}.Log.final.out"), emit: read_stats

    script:
    """
    STAR --genomeDir ${star_index} \\
        --outFileNamePrefix ${name}. \\
        --readFilesCommand zcat \\
        --readFilesIn ${fastq1} ${fastq2} \\
        --outFilterMultimapNmax 1000 \\
        --outSAMmultNmax 1 \\
        --chimSegmentMin 10 \\
        --chimJunctionOverhangMin 10 \\
        --chimOutJunctionFormat 1 \\
        --alignSJDBoverhangMin 10 \\
        --alignMatesGapMax 200000 \\
        --alignIntronMax 200000 \\
        --chimSegmentReadGapMax 3 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --seedSearchStartLmax 20 \\
        --winAnchorMultimapNmax 50 \\
        --outSAMtype BAM Unsorted \\
        --chimOutType Junctions WithinBAM \\
        --outSAMunmapped Within KeepPairs \\
        --runThreadN ${task.cpus}

    mv ${name}.Aligned.out.bam ${name}.bam
    """
}


process STAR_ARRIBA {
    tag "${name}"
    label 'process_medium'

    conda ("${baseDir}/environments/alignment.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9b/9b8ecb2f9a77b5e7573ef6fae2f4c2e771064f7a129ed1329913c1025c33f365/data' }"

    input:
      tuple val(name), path(fastq1), path(fastq2)
      path(star_index, stageAs: "star_index/")

    output:
      tuple val("${name}"), path("${name}.bam"), emit: bams

    script:
    """
    STAR --genomeDir "${star_index}" \\
        --outFileNamePrefix ${name}. \\
        --readFilesCommand zcat \\
        --readFilesIn ${fastq1} ${fastq2} \\
        --outFilterMultimapNmax 50 \\
        --peOverlapNbasesMin 10 \\
        --alignSplicedMateMapLminOverLmate 0.5 \\
        --chimSegmentMin 10 \\
        --chimJunctionOverhangMin 10 \\
        --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreSeparation 1 \\
        --chimSegmentReadGapMax 3 \\
        --chimMultimapNmax 50 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --outSAMtype BAM Unsorted \\
        --chimOutType WithinBAM HardClip \\
        --outSAMunmapped Within \\
        --runThreadN ${task.cpus}

    mv ${name}.Aligned.out.bam ${name}.bam
    """
}


process READ_FILTER {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/filtering.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pysam:0.22.0--a94c5bab35035aad' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64682c99fc92227f78f81a53c8d739b16e8b712c6c75a8909f159405cb29dbe1/data' }"

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
    label 'process_low'

    conda ("${baseDir}/environments/samtools.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/samtools:1.9--ff876e25d460de68' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/597f0179350adad0f4fa6b045f64f469fc04bf6c1a26c870bf900f2832bf260f/data' }"

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