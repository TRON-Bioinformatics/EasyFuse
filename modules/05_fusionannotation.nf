process FUSION_ANNOTATION {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/annotation.yml")

    input:
      tuple val(name), path(fusions)

    output:
      tuple val("${name}"), path("annotated_fusions.csv"), path("annotated_fusions.csv.debug"), path("annotated_fusions.csv.fasta"), emit: annot_fusions
      
      

    script:

    """
    fusionannotator.py \
        --detected_fusions ${fusions} \
        --annotation_db ${params.annotation_db} \
        --out_csv annotated_fusions.csv \
        --genome_fasta ${params.fasta} \
        --tsl_info ${params.reference_tsl} \
        --cis_near_dist 1000000 \
        --context_seq_len 400 \
        --tsl_filter_level 4,5,NA
    """
}
