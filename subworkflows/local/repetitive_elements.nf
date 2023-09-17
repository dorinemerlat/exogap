include {REPEATMODELER2                                             } from '../../modules/local/repeatmodeler2/repeatmodeler2'
include {RM2_REFORMAT                                               } from '../../modules/local/repeatmodeler2/rm2_reformat'
include {LIBRARY_CONCATENATION                                      } from '../../modules/local/repeatmodeler2/library_concatenation'
include {LIBRARY_SPLIT                                              } from '../../modules/local/repeatmodeler2/library_split'
include {ELIMINATE_REDUNDANCE   as ELIMINATE_REDUNDANCE_1           } from '../../modules/local/repeatmodeler2/eliminate_redundance'
include {ELIMINATE_REDUNDANCE   as ELIMINATE_REDUNDANCE_2           } from '../../modules/local/repeatmodeler2/eliminate_redundance'
include {REPEATMASKER           as REPEATMASKER_WITH_EXISTING_LIB   } from '../../modules/local/repeatmasker/repeatmasker'
include {REPEATMASKER           as REPEATMASKER_WITH_OWN_LIB_1      } from '../../modules/local/repeatmasker/repeatmasker'
include {REPEATMASKER           as REPEATMASKER_WITH_OWN_LIB_2      } from '../../modules/local/repeatmasker/repeatmasker'
// include {PROCESS_REPEATS                                            } from '../../modules/local/repeatmasker/process_repeats'
// include {RM_REFORMAT                                                } from '../../modules/local/repeatmasker/rm_output_reformat'
// include {PASTEC                                                     } from '../../modules/local/pastec'
// include {STATS                                                      } from '../../modules/local/re_stats/stats'
// include {PLOTS                                                      } from '../../modules/local/re_stats/plots'

workflow REPETITIVE_ELEMENTS{
    take:
        ch_genomes
        ch_repeats_lib

    main:
        // TO DO: reuse as principal instruction
        // REPEATMODELER2(ch_genomes)

        // ch_RM2_library = RM2_REFORMAT(REPEATMODELER2.out)

        // TO Do: remove
        ch_RM2_library = Channel.fromPath("/gstock/user/merlat/myriapods/repetitive_elements/*.fa", checkIfExists: true)
            .map { file -> tuple(file.baseName, file) }

        ch_RM2_library = RM2_REFORMAT(ch_RM2_library)

        // concatenate libraries == yes
        if ( params.repeats_concat ) {
            ch_RM2_library = LIBRARY_CONCATENATION(
                    ch_RM2_library.map{ it -> it.last() }.collect())
            }

        // split library == yes
        if ( params.repeats_split ) {
            LIBRARY_SPLIT(ch_RM2_library)
            ELIMINATE_REDUNDANCE_1(LIBRARY_SPLIT.out.classified)
            ELIMINATE_REDUNDANCE_2(LIBRARY_SPLIT.out.unclassified)
        }
        else {
            ELIMINATE_REDUNDANCE_1(ch_RM2_library)
        }


        // use a extern repeats library == yes
        ch_genomes.view()
        if (params.repeats_lib ){
            ch_repeatmasker = REPEATMASKER_WITH_EXISTING_LIB(ch_genomes, ch_repeats_lib)
            ch_masked = ch_repeatmasker.masked
            ch_cat = ch_repeatmasker.cat
        }
        // use a extern repeats library == no
        else {
            ch_masked = ch_genomes
            ch_cat = Channel.of('')
        }

        // in all case: make one iteration with own library (all or specie-specific)
        ch_repeatmasker = REPEATMASKER_WITH_OWN_LIB_1(
                ch_masked,
                ELIMINATE_REDUNDANCE_1.out)

        ch_cat = ch_cat.concat(ch_repeatmasker.cat)

        // split repeats == yes
        if ( params.repeats_split ) {
            ch_repeatmasker = REPEATMASKER_WITH_OWN_LIB_2(
                    ch_repeatmasker.masked,
                    ELIMINATE_REDUNDANCE_2.out)

            ch_cat = ch_cat.concat(ch_repeatmasker.cat)
        }

        ch_cat = ch_cat.collect()

        TO DO: add PROCESS_REPEATS process
            PROCESS_REPEATS(
                    ch_repeatmasker.masked,
                    ch_repeatmasker.out,
                    ch_cat,
                    ch_repeats_lib

            )

        // // TO DO: add RM_OUTPUT_REFORMAT process
        //     RM_REFORMAT(PROCESS_REPEATS.out)

        // // TO DO: add PASTEC process
        //     PASTEC(RM_REFORMAT.out)

        // // TO DO: add RE_STATISTICS process
        //     STATS(PASTEC.out)
        //     PLOTS(PASTEC.out)

    // emit:
    //     re = ELIMINATE_REDUNDANCE_1.output
        // re_stats = STATS.out
        // re_gff3 = repeats_ch.gff3
        // masked_genomes = repeats_ch.fasta
}

