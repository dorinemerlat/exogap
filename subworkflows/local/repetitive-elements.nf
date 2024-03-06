include {REPEATMODELER                                              } from '../../modules/local/repeatmodeler/repeatmodeler'
include {RENAME_REPEATMODELER_OUTPUT                                } from '../../modules/local/repeatmodeler/rename-repeatmodeler-output'
include {GATHER_FILES                                               } from '../../modules/local/gather-files'
include {SEPARATE_LIBRARIES                                         } from '../../modules/local/bioawk/separate-libraries'
include {CD_HIT_FOR_REPEATS as CD_HIT_FOR_REPEATS_1                 } from '../../modules/local/cd-hit/cd-hit-for-repeats'
include {CD_HIT_FOR_REPEATS as CD_HIT_FOR_REPEATS_2                 } from '../../modules/local/cd-hit/cd-hit-for-repeats'
include {REPEATMASKER           as REPEATMASKER_WITH_EXISTING_LIB   } from '../../modules/local/repeatmasker/repeatmasker'
include {REPEATMASKER           as REPEATMASKER_WITH_OWN_LIB_1      } from '../../modules/local/repeatmasker/repeatmasker'
include {REPEATMASKER           as REPEATMASKER_WITH_OWN_LIB_2      } from '../../modules/local/repeatmasker/repeatmasker'
include {PROCESS_REPEATS                                            } from '../../modules/local/repeatmasker/processrepeats'
include {REPEATLANDSCAPE                                            } from '../../modules/local/repeatmasker/repeatLanscape'
include {SUMMARIZE_REPEATS                                          } from '../../modules/local/stats-and-plots/summarize-repeats'
include {PLOT_REPEATS                                               } from '../../modules/local/stats-and-plots/plot_repeats.nf'
include {GET_DFAM_CLASSIFICATION                                    } from '../../modules/local/api/get_dfam_classification.nf'

workflow REPETITIVE_ELEMENTS {
    take:
        genomes
        newick

    main:
        // Run RepeatModeler and process results
        REPEATMODELER(genomes)
        RENAME_REPEATMODELER_OUTPUT(REPEATMODELER.out)
            .set { rm_library }

        // use one library for all genomes
        if ( params.one_library_for_all ) {
            GATHER_FILES(rm_library.map{ it -> tuple('families.fa', it.last()) }.groupTuple())
                .set { rm_library }
            }

        SEPARATE_LIBRARIES(rm_library)

        CD_HIT_FOR_REPEATS_1(SEPARATE_LIBRARIES.out.classified)
        CD_HIT_FOR_REPEATS_2(SEPARATE_LIBRARIES.out.unclassified)

        // use a extern repeats library
        Channel.fromPath( params.external_library, checkIfExists: true )
            .combine(genomes.map {meta, file -> meta })
            .map { library, meta -> [ meta, library ] }
            .set { external_library }

        REPEATMASKER_WITH_EXISTING_LIB(genomes, external_library)

        // make 2 iterations with own library (all or specie-specific)
        REPEATMASKER_WITH_OWN_LIB_1(REPEATMASKER_WITH_EXISTING_LIB.out.masked, CD_HIT_FOR_REPEATS_1.out)
        REPEATMASKER_WITH_OWN_LIB_2(REPEATMASKER_WITH_OWN_LIB_1.out.masked, CD_HIT_FOR_REPEATS_2.out)

        // Group all .cat files
        REPEATMASKER_WITH_EXISTING_LIB.out.cat
            .concat(REPEATMASKER_WITH_OWN_LIB_1.out.cat, REPEATMASKER_WITH_OWN_LIB_2.out.cat)
            .groupTuple()
            .set { cat }

        // Process repeats and landscape analysis
        PROCESSREPEATS(REPEATMASKER_WITH_OWN_LIB_2.masked, cat, external_library)


        REPEATLANDSCAPE(PROCESSREPEATS.out.align)

        GET_DFAM_CLASSIFICATION()
        CHANGE_CLASSIFICATION_TO_DFAM(
            GET_DFAM_CLASSIFICATION,
            PROCESSREPEATS.out.masked,
            PROCESSREPEATS.out.gff,
            PROCESSREPEATS.out.cat,
            PROCESSREPEATS.out.out,
            PROCESSREPEATS.out.tbl,
            REPEATLANDSCAPE.out.align)

        GET_COMPLEX_REPEATS(CHANGE_CLASSIFICATION_TO_DFAM.out.gff)
        SUMMARIZE_AND_PLOT_REPEATS(CHANGE_CLASSIFICATION_TO_DFAM.out.out.collect{it[1]}, newick)

    // emit:
    // re_gff3 = repeats_ch.gff3
    // masked_genomes = repeats_ch.fasta
}
