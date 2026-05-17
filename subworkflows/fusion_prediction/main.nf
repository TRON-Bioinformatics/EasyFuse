//
// EasyFuse: Fusion Detection - subworkflow to detect gene fusions in RNA-seq data using tools of choice
//
//

include { ARRIBA        } from '../../modules/arriba/main'
include { STARFUSION    } from '../../modules/starfusion/main'
include { STAR_ARRIBA   } from '../../modules/stararriba/main'
include { FUSIONCATCHER } from '../../modules/fusioncatcher/main'

include { PARSE_ARRIBA        } from '../../modules/fusionparsing/parsearriba/main'
include { PARSE_STARFUSION    } from '../../modules/fusionparsing/parsestarfusion/main'
include { PARSE_FUSIONCATCHER } from '../../modules/fusionparsing/parsefusioncatcher/main'

workflow FUSION_PREDICTION {

    take:
    ch_filtered_fastqs     // channel: [ val(meta), [ fastqs ] ]
    // ch_chimeric_reads      // channel: [ val(meta), [ chimeric_reads ] ]
    ch_fusioncatcher_index // channel: [ path(fusioncatcher_index) ]
    ch_starfusion_index    // channel: [ path(starfusion_index)]
    ch_stararriba_index    // channel: [ path(stararriba_index) ]
    ch_reference_gtf       // channel: [ path(reference_gtf) ]
    ch_reference_fasta     // channel: [ path(reference_fasta) ]
    ch_fusiontools         // channel: [ [arriba: true, fusioncatcher: true, starfusion: true] ]

    main:

    ch_versions = channel.empty()
    ch_arriba_fusions = channel.empty()
    ch_starfusion_fusions = channel.empty()
    ch_fusioncatcher_fusions = channel.empty()

    // run fusionctacher and generate parsed output
    if (ch_fusiontools.value.run_fusioncatcher) {

        ch_fusioncatcher_input = ch_filtered_fastqs.combine(ch_fusioncatcher_index)
        FUSIONCATCHER(ch_fusioncatcher_input)
        ch_versions = ch_versions.mix(FUSIONCATCHER.out.versions)

        PARSE_FUSIONCATCHER(FUSIONCATCHER.out.fusions)
        ch_versions = ch_versions.mix(PARSE_FUSIONCATCHER.out.versions)
        ch_fusioncatcher_fusions = PARSE_FUSIONCATCHER.out.fusions
    }

    // run starfusion and generate parsed output
    if (ch_fusiontools.value.run_starfusion) {

        ch_starfusion_input = ch_filtered_fastqs.combine(ch_starfusion_index)
        STARFUSION(ch_starfusion_input)
        ch_versions = ch_versions.mix(STARFUSION.out.versions)

        PARSE_STARFUSION(STARFUSION.out.fusions)
        ch_versions = ch_versions.mix(PARSE_STARFUSION.out.versions)
        ch_starfusion_fusions = PARSE_STARFUSION.out.fusions
    }

    // run stararriba and generate parsed output
    if (ch_fusiontools.value.run_arriba) {

        ch_stararriba_input = ch_filtered_fastqs.combine(ch_stararriba_index)
        STAR_ARRIBA(ch_stararriba_input)
        ch_versions = ch_versions.mix(STAR_ARRIBA.out.versions)

        ch_arriba_input = STAR_ARRIBA.out.bam
                            .combine(ch_reference_gtf)
                            .combine(ch_reference_fasta)
                            .map { meta, bam, ref_gtf, ref_fasta -> tuple(meta, bam, ref_gtf, ref_fasta) }
        ARRIBA(ch_arriba_input)
        ch_versions = ch_versions.mix(ARRIBA.out.versions)

        PARSE_ARRIBA(ARRIBA.out.fusions)
        ch_versions = ch_versions.mix(PARSE_ARRIBA.out.versions)
        ch_arriba_fusions = PARSE_ARRIBA.out.fusions
    }

    emit:

    arriba_results = ch_arriba_fusions               // channel: [ tuple(val(meta), path(arriba_fusions)) ]
    starfusion_results = ch_starfusion_fusions       // channel: [ tuple(val(meta), path(starfusion_fusions)) ]
    fusioncatcher_results = ch_fusioncatcher_fusions // channel: [ tuple(val(meta), path(fusioncatcher_fusions)) ]

    versions = ch_versions                           // channel: [ versions.yml ]
}
