
process FUSION_PARSER {
    cpus 1
    memory "8g"
    tag "${name}"
    //publishDir "${params.output}/${name}", mode: 'copy'

    conda ("${baseDir}/environments/filtering.yml")

    input:
      tuple val(name), path(fusion_catcher_1), path(fusion_catcher_2), path(star_fusion), path(arriba_1), path(arriba_2)


    output:
      tuple val("${name}"), path("Detected_Fusions.csv"), emit: fusions

    script:

    """
    fusiontoolparser.py \
        --input_fusioncatcher ${fusion_catcher_1} \
        --input_fusioncatcher2 ${fusion_catcher_2} \
	      --input_starfusion ${star_fusion} \
        --input_arriba ${arriba_1} \
        --output . \
        --sample ${name}
    """
}

process FUSION_ANNOTATION {
    cpus 1
    memory "8g"
    tag "${name}"
    //publishDir "${params.output}/${name}", mode: 'copy'

    conda ("${baseDir}/environments/annotation.yml")

    input:
      tuple val(name), path(fusions)

    output:
      tuple val("${name}"), path("annotated_fusions.csv"), path("annotated_fusions.csv.debug"), path("annotated_fusions.csv.fasta"), emit: annot_fusions
      
      

    script:

    """
    fusionannotation.py \
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
