include {REPEATMODELER                                              } from '../../modules/local/repeatmodeler/repeatmodeler'
include {RENAME_REPEATMODELER_OUTPUT                                } from '../../modules/local/repeatmodeler/rename-repeatmodeler-output'
include {GATHER_LIBRARIES                                           } from '../../modules/local/repeatmodeler/gather-libraries'
include {SEPARATE_LIBRARIES                                         } from '../../modules/local/bioawk/separate-libraries'
include {CD_HIT_FOR_REPEATS as CD_HIT_FOR_REPEATS_1                 } from '../../modules/local/cd-hit/cd-hit-for-repeats'
include {CD_HIT_FOR_REPEATS as CD_HIT_FOR_REPEATS_2                 } from '../../modules/local/cd-hit/cd-hit-for-repeats'
include {REPEATMASKER           as REPEATMASKER_WITH_EXISTING_LIB   } from '../../modules/local/repeatmasker/repeatmasker'
include {REPEATMASKER           as REPEATMASKER_WITH_OWN_LIB_1      } from '../../modules/local/repeatmasker/repeatmasker'
include {REPEATMASKER           as REPEATMASKER_WITH_OWN_LIB_2      } from '../../modules/local/repeatmasker/repeatmasker'
// include {PROCESS_REPEATS                                            } from '../../modules/local/repeatmasker/process-repeats'
// include {PASTEC                                                     } from '../../modules/local/pastec'
// include {STATS                                                      } from '../../modules/local/re_stats/stats'
// include {PLOTS                                                      } from '../../modules/local/re_stats/plots'

workflow ANNOTATE_REPEATS {
    take:
        genomes
        repeats_lib

    main:
        // TO DO: reuse as principal instruction
        // REPEATMODELER(genomes)

        // rm_library = RENAME_REPEATMODELER_OUTPUT(REPEATMODELER.out)

        // TO Do: remove
        rm_library = Channel.fromPath("/gstock/user/merlat/myriapods/repetitive_elements/*.fa", checkIfExists: true)
            .map { file -> tuple(file.baseName, file) }

        rm_library = RENAME_REPEATMODELER_OUTPUT(rm_library)

        // concatenate libraries == yes
        if ( params.repeats_concat ) {
            rm_library = GATHER_LIBRARIES(
                    rm_library.map{ it -> it.last() }.collect())
            }

        // split library == yes
        if ( params.repeats_split ) {
            SEPARATE_LIBRARIES(rm_library)
            CD_HIT_FOR_REPEATS_1(SEPARATE_LIBRARIES.out.classified)
            CD_HIT_FOR_REPEATS_2(SEPARATE_LIBRARIES.out.unclassified)
        }
        else {
            CD_HIT_FOR_REPEATS_1(rm_library)
        }

        // use a extern repeats library == yes
        if (params.repeats_lib ){
            repeatmasker = REPEATMASKER_WITH_EXISTING_LIB(genomes, repeats_lib)
            masked = repeatmasker.masked
            cat = repeatmasker.cat
        }
        // use a extern repeats library == no
        else {
            masked = genomes
            cat = Channel.of('')
        }

        // in all case: make one iteration with own library (all or specie-specific)
        repeatmasker = REPEATMASKER_WITH_OWN_LIB_1(
                masked.last(),
                CD_HIT_FOR_REPEATS_1.out)

        cat = cat.concat(repeatmasker.cat)

        // split repeats == yes
        if ( params.repeats_split ) {
            repeatmasker = REPEATMASKER_WITH_OWN_LIB_2(
                    REPEATMASKER_WITH_OWN_LIB_1.out.masked,
                    CD_HIT_FOR_REPEATS_2.out)

            cat = cat.concat(repeatmasker.cat)
        }

        cat = cat.collect()

        // TO DO: add PROCESS_REPEATS process
        // PROCESS_REPEATS(
        //         repeatmasker.masked,
        //         repeatmasker.out,
        //         cat,
        //         repeats_lib
        // )

        // // TO DO: add RM_OUTPUT_REFORMAT process
        //     RM_REFORMAT(PROCESS_REPEATS.out)

        // // TO DO: add PASTEC process
        //     PASTEC(RM_REFORMAT.out)

        // // TO DO: add RE_STATISTICS process
        //     STATS(PASTEC.out)
        //     PLOTS(PASTEC.out)

    // emit:
    //     re = CD_HIT_FOR_REPEATS_1.output
        // re_stats = STATS.out
        // re_gff3 = repeats_ch.gff3
        // masked_genomes = repeats_ch.fasta
}
