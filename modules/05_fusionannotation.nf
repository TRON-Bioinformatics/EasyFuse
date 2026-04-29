process FUSION_ANNOTATION {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/annotation.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/gffutils_biopython_python-xxhash_python:458dcfe93321068c' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/30840bab8184dfe45a8e2d0f41c869e5fcc447c4e83f0c9db5e74e9e3c07fa87/data' }"

    input:
      tuple val(name), path(fusions)
      path(annotation_db)

    output:
      tuple val("${name}"), path("annotated_fusions.csv"), path("annotated_fusions.csv.debug"), path("annotated_fusions.csv.fasta"), emit: annot_fusions

    script:
    """
    fusionannotator.py \
        --detected_fusions ${fusions} \
        --annotation_db ${annotation_db} \
        --out_csv annotated_fusions.csv \
        --genome_fasta ${params.fasta} \
        --tsl_info ${params.reference_tsl} \
        --cis_near_dist 1000000 \
        --context_seq_len 400 \
        --tsl_filter_level 4,5,NA
    """
}
