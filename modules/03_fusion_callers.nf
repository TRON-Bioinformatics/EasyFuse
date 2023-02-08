

process MAPSPLICE {
    cpus 6
    memory "32g"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::mapsplice=2.2.1" : null)

    input:
      tuple val(name), path(fastq1), file(fastq2)

    output:
      tuple val("${name}"), path("fusions_candidates.txt"), emit: fusions

    script:
    """
    mapsplice.py \
    --chromosome-dir ${params.chromsome_dir} \
    -x ${params.bowtie_index} \
    -1 <(gunzip -c -f ${fastq1}) \
    -2 <(gunzip -c -f ${fastq2}) \
    --threads ${task.cpus} \
    --output . \
    --qual-scale phred33 \
    --bam \
    --seglen 20 \
    --min-map-len 40 \
    --gene-gtf ${params.gtf}

    mv ${name}.Aligned.out.bam ${name}.bam
    """
}