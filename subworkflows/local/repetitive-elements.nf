include {REPEATMODELER                                              } from '../../modules/local/repeatmodeler/repeatmodeler'
include {RENAME_REPEATMODELER_OUTPUT                                } from '../../modules/local/repeatmodeler/rename-repeatmodeler-output'
include {GATHER_FILES                                               } from '../../modules/local/gather-files'
include {SEPARATE_LIBRARIES                                         } from '../../modules/local/bioawk/separate-libraries'
include {CD_HIT_FOR_REPEATS as CD_HIT_FOR_REPEATS_1                 } from '../../modules/local/cd-hit/cd-hit-for-repeats'
include {CD_HIT_FOR_REPEATS as CD_HIT_FOR_REPEATS_2                 } from '../../modules/local/cd-hit/cd-hit-for-repeats'
include {REPEATMASKER           as REPEATMASKER_WITH_EXISTING_LIB   } from '../../modules/local/repeatmasker/repeatmasker'
include {REPEATMASKER           as REPEATMASKER_WITH_OWN_LIB_1      } from '../../modules/local/repeatmasker/repeatmasker'
include {REPEATMASKER           as REPEATMASKER_WITH_OWN_LIB_2      } from '../../modules/local/repeatmasker/repeatmasker'
include {PROCESS_REPEATS                                            } from '../../modules/local/repeatmasker/process-repeats'
include {REPEATLANDSCAPE                                            } from '../../modules/local/repeatmasker/repeatLanscape'
// include {STATS                                                      } from '../../modules/local/re_stats/stats'
// include {PLOTS                                                      } from '../../modules/local/re_stats/plots'

workflow REPETITIVE_ELEMENTS {
    take:
        genomes

    main:
        // // Run RepeatModeler and process results
        // REPEATMODELER(genomes)

        // RENAME_REPEATMODELER_OUTPUT(REPEATMODELER.out)
        //     .set { rm_library }

        // // use one library for all genomes
        // if ( params.one_library_for_all ) {
        //     GATHER_FILES(rm_library.map{ it -> tuple('families.fa', it.last()) }.groupTuple())
        //         .set { rm_library }
        //     }

        // SEPARATE_LIBRARIES(rm_library)

        // CD_HIT_FOR_REPEATS_1(SEPARATE_LIBRARIES.out.classified)
        // CD_HIT_FOR_REPEATS_2(SEPARATE_LIBRARIES.out.unclassified)

        // // use a extern repeats library
        // Channel.fromPath( params.external_library, checkIfExists: true )
        //     .combine(genomes.map {meta, file -> meta })
        //     .map { library, meta -> [ meta, library ] }
        //     .set { external_library }

        // REPEATMASKER_WITH_EXISTING_LIB(genomes, external_library)

        // REPEATMASKER_WITH_EXISTING_LIB.out.masked.set { masked }
        // REPEATMASKER_WITH_EXISTING_LIB.out.cat.set { cat }

        // // dont't use a extern repeats library
        // else {
        //     genomes.set { masked }
        //     genomes.map { meta, fasta -> [ meta, '' ] }.set { cat }
        // }

        // // in all case: make one iteration with own library (all or specie-specific)
        // REPEATMASKER_WITH_OWN_LIB_1(masked, CD_HIT_FOR_REPEATS_1.out)
        // REPEATMASKER_WITH_OWN_LIB_1.out.cat.join(cat).set { cat }

        // REPEATMASKER_WITH_OWN_LIB_2(REPEATMASKER_WITH_OWN_LIB_1.out.masked, CD_HIT_FOR_REPEATS_2.out)
        // REPEATMASKER_WITH_OWN_LIB_2.out.cat.join(cat).set { cat }

        Channel.fromPath("/gstock/user/merlat/myriapods/repetitive_elements/repeatmasker1/*cat")
            .map { file -> tuple([id : file.simpleName], file) }
            .set { cat1}

        Channel.fromPath("/gstock/user/merlat/myriapods/repetitive_elements/repeatmasker2/*cat")
            .map { file -> tuple([id : file.simpleName], file) }
            .set { cat2}

        Channel.fromPath("/gstock/user/merlat/myriapods/repetitive_elements/repeatmasker3/*cat")
            .map { file -> tuple([id : file.simpleName], file) }
            .set { cat3}

        cat1.concat(cat2,cat3).groupTuple().set {cat_final}

        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .join( cat_final.map{ meta, files -> [ meta.id, files ] } )
            .map { id, meta, genome, files -> [ meta, files ] }
            .set {cat_final}

        Channel.fromPath("/gstock/user/merlat/myriapods/repetitive_elements/repeatmasker3/*masked")
            .map { file -> tuple([id : file.simpleName], file) }
            .set { masked_final }

        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .join( masked_final.map{ meta, file -> [ meta.id, file ] } )
            .map { id, meta, genome, file -> [ meta, file ] }
            .set { masked_final }

        Channel.fromPath( params.external_library, checkIfExists: true )
            .combine(masked_final.map {meta, file -> meta })
            .map { library, meta -> [ meta, library ] }
            .set { external_library }

        // Process repeats and landscape analysis
        PROCESS_REPEATS(masked_final, cat_final, external_library )

        REPEATLANDSCAPE(PROCESS_REPEATS.out.align)

        // // Process repeats and landscape analysis
        // PROCESS_REPEATS(REPEATMASKER_WITH_OWN_LIB_2.masked, cat, external_library)
        // REPEATLANDSCAPE(PROCESS_REPEATS.out.align)

        // PLOT_REPEATLANDSCAPE(REPEATLANDSCAPE.out) -> one by genome
        // PlOT_REPEAT_OVERVIEW() -> one for all
        // PLOT_CORRELATION() -> one for each group
        // PLOT_ABUNDANCE() -> one for all

    // emit:
    // re_gff3 = repeats_ch.gff3
    // masked_genomes = repeats_ch.fasta
}
