params.reference = false
params.chromosome_dir = "/projects/data/human/ensembl/GRCh38.86/fasta"
params.bowtie_index = "/projects/data/human/ensembl/GRCh38.86/bowtie_index/hg38"
params.gtf = "/projects/data/human/ensembl/GRCh38.86/Homo_sapiens.GRCh38.86.gtf"
params.fusioncatcher_index = "/scratch/info/data/easyfuse/easyfuse_ref/fusioncatcher_index/"
params.starfusion_index = "/projects/data/human/ensembl/GRCh38.86/starfusion_index/"

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
    gunzip -c -f ${fastq1} > ${name}_1.fastq
    gunzip -c -f ${fastq2} > ${name}_2.fastq

    mapsplice.py \
    --chromosome-dir ${params.chromosome_dir} \
    -x ${params.bowtie_index} \
    -1 ${name}_1.fastq \
    -2 ${name}_2.fastq \
    --threads ${task.cpus} \
    --output . \
    --qual-scale phred33 \
    --bam \
    --seglen 20 \
    --min-map-len 40 \
    --gene-gtf ${params.gtf}

    rm -f ${name}_1.fastq
    rm -f ${name}_2.fastq
    """
}

process FUSION_CATCHER {
    cpus 6
    memory "32g"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::fusioncatcher=1.33" : null)

    input:
      tuple val(name), path(fastq1), file(fastq2)

    output:
      tuple val("${name}"), path("final-list_candidate-fusion-genes.txt"), emit: fusions

    script:
    """
    fusioncatcher \
        --input ${fastq1},${fastq2} \
        --data ${params.fusioncatcher_index} \
        --output . \
        -p ${task.cpus}
    """
}

process STAR_FUSION {
    cpus 6
    memory "32g"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::star-fusion=1.12.0" : null)

    input:
      tuple val(name), path(chimeric_reads)

    output:
      tuple val("${name}"), path("star-fusion.fusion_predictions.tsv"), emit: fusions

    script:
    """
    STAR-Fusion \
        --chimeric_junction ${chimeric_reads} \
        --genome_lib_dir ${params.starfusion_index} \
        --CPU ${task.cpus} \
        --output_dir .
    """
}