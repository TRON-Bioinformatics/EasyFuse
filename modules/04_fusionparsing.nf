process PARSE_ARRIBA {
    tag "${name}"
    label 'process_single'

    conda ("${baseDir}/environments/fusionparsing.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python:3.8.0--5e0e57f6a223cdda' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/796c2d0f22be8bd23bf91dd74661318db5545bc15c4ec63fd1207c95b2d5d22c/data' }" 

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

    conda ("${baseDir}/environments/fusionparsing.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python:3.8.0--5e0e57f6a223cdda' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/796c2d0f22be8bd23bf91dd74661318db5545bc15c4ec63fd1207c95b2d5d22c/data' }"

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

    conda ("${baseDir}/environments/fusionparsing.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python:3.8.0--5e0e57f6a223cdda' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/796c2d0f22be8bd23bf91dd74661318db5545bc15c4ec63fd1207c95b2d5d22c/data' }"

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

    conda ("${baseDir}/environments/fusionparsing.yml")
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python:3.8.0--5e0e57f6a223cdda' :
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/796c2d0f22be8bd23bf91dd74661318db5545bc15c4ec63fd1207c95b2d5d22c/data' }"

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