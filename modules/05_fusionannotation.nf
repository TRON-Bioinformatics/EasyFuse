
process PARSE_ARRIBA {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/filtering.yml")

    input:
      tuple val(name), path(arriba_out)

    output:
      tuple val("${name}"), path("arriba.csv"), emit: fusions

    script:

    """
    parse_tool.py \
        --input_file ${arriba_out} \
        --output_file arriba.csv \
        --tool arriba
    """
}

process PARSE_STAR_FUSION {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/filtering.yml")

    input:
      tuple val(name), path(star_fusion_out)

    output:
      tuple val("${name}"), path("starfusion.csv"), emit: fusions

    script:

    """
    parse_tool.py \
        --input_file ${star_fusion_out} \
        --output_file starfusion.csv \
        --tool starfusion
    """
}

process PARSE_FUSION_CATCHER {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/filtering.yml")

    input:
      tuple val(name), path(fusion_catcher_1), path(fusion_catcher_2)

    output:
      tuple val("${name}"), path("fusioncatcher.csv"), emit: fusions

    script:

    """
    parse_tool.py \
        --input_file ${fusion_catcher_1} \
        --input_file2 ${fusion_catcher_2} \
        --output_file fusioncatcher.csv \
        --tool fusioncatcher
    """
}

process FUSION_PARSER {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/filtering.yml")

    input:
      tuple val(name), path(fusion_catcher), path(star_fusion), path(arriba)


    output:
      tuple val("${name}"), path("Detected_Fusions.csv"), emit: fusions

    script:

    """
    fusiontoolparser.py \
        --tool fusioncatcher ${fusion_catcher} \
	      --tool starfusion ${star_fusion} \
        --tool arriba ${arriba} \
        --output Detected_Fusions.csv \
        --sample ${name}
    """
}

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
