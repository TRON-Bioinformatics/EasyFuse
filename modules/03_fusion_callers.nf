params.reference = false
params.chromosome_dir = "/projects/data/human/ensembl/GRCh38.86/fasta"
params.bowtie_index = "/projects/data/human/ensembl/GRCh38.86/bowtie_index/hg38"
params.gtf = "/projects/data/human/ensembl/GRCh38.86/Homo_sapiens.GRCh38.86.gtf"
params.fusioncatcher_index = "/scratch/info/data/easyfuse/easyfuse_ref/fusioncatcher_index/"
params.starfusion_index = "/projects/data/human/ensembl/GRCh38.86/starfusion_index/"


process FUSION_CATCHER_INDEX {
    cpus 6
    memory "32g"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::fusioncatcher=1.33" : null)

    script:
    """
    download-human-db.sh
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
      tuple val("${name}"), path("summary_candidate_fusions.txt"), path("final-list_candidate-fusion-genes.txt"), emit: fusions

    script:
    """
    # --data ${params.fusioncatcher_index} \
    fusioncatcher \
        --input ${fastq1},${fastq2} \
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