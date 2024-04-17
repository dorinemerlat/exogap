include { REPEATMODELER                                     } from '../../modules/local/repeatmodeler'
include { RENAME_REPEATMODELER_OUTPUT                       } from '../../modules/local/rename_repeatmodeler_output'
include { GATHER_FILES                                      } from '../../modules/local/gather_files'
include { SEPARATE_LIBRARIES                                } from '../../modules/local/separate_libraries'
include { CD_HIT_EST as CD_HIT_EST_1                        } from '../../modules/local/cd_hit_est'
include { CD_HIT_EST as CD_HIT_EST_2                        } from '../../modules/local/cd_hit_est'
include { REPEATMASKER as REPEATMASKER_WITH_EXISTING_LIB    } from '../../modules/local/repeatmasker'
include { REPEATMASKER as REPEATMASKER_WITH_OWN_LIB_1       }from '../../modules/local/repeatmasker'
include { REPEATMASKER as REPEATMASKER_WITH_OWN_LIB_2       } from '../../modules/local/repeatmasker'
include { PROCESS_REPEATS                                   } from '../../modules/local/processrepeats'
include { REPEATLANDSCAPE                                   } from '../../modules/local/repeatLanscape'
include { SUMMARIZE_REPEATS                                 } from '../../modules/local/summarize_repeats'
// include {PLOT_REPEATS }                                     from '../../modules/local/plot_repeats.nf'
include { DOWNLOAD_DFAM                                     } from '../../modules/local/download_dfam.nf'
include { REFORMAT_CLASSIFICATION_TO_DFAM                   } from '../../modules/local/reformat_classification_to_dfam.nf'

workflow ANNOTATE_REPETITIVE_ELEMENTS {
    take:
        genomes
        // newick

    main:
        // Run RepeatModeler and process results
        genomes.map { id, meta, fasta -> [ id, meta + ['genome': fasta], fasta ] }
            .set { genomes }

        REPEATMODELER(genomes)
        RENAME_REPEATMODELER_OUTPUT(REPEATMODELER.out)
            .set { rm_library }

        // use one library for all genomes
        if ( params.one_library_for_all ) {
            GATHER_FILES(rm_library.map{ id, meta, file -> [ 'library_repeats', id, file ] }
                                    .groupTuple()
                                    .map {id, meta, file -> [ id, ['species': meta, 'name': 'repeatmodeler'], file, 'library_repeats.fa', 'no' ]})
                .set { rm_library }
            }

        SEPARATE_LIBRARIES(rm_library)

        CD_HIT_EST_1(SEPARATE_LIBRARIES.out.classified.map { id, meta, file -> [ id, meta, file, '10', '80', '98' ] })
        CD_HIT_EST_2(SEPARATE_LIBRARIES.out.unclassified.map { id, meta, file -> [ id, meta, file, '10', '80', '98' ] })

        // use a extern repeats library
        if ( params.external_library ) {
            REPEATMASKER_WITH_EXISTING_LIB(genomes, Channel.fromPath( params.external_library, checkIfExists: true ) )
            REPEATMASKER_WITH_EXISTING_LIB.out.masked
                .set { masked_genomes }
        } else {
            masked_genomes = genomes
        }

        // make 2 iterations with own library (all or specie-specific)
        REPEATMASKER_WITH_OWN_LIB_1(masked_genomes, CD_HIT_EST_1.out.map { id, meta, file -> file } )
        REPEATMASKER_WITH_OWN_LIB_2(REPEATMASKER_WITH_OWN_LIB_1.out.masked, CD_HIT_EST_2.out.map { id, meta, file -> file } )

        // Group all .cat files
        REPEATMASKER_WITH_OWN_LIB_1.out.cat
            .combine(REPEATMASKER_WITH_OWN_LIB_2.out.cat, by: 0)
            .map { id, meta1, cat1, meta2, cat2 -> [ id, meta1, [cat1, cat2] ]}
            .set { all_cat }

        if ( params.external_library ) {
            all_cat.combine(REPEATMASKER_WITH_EXISTING_LIB.out.cat, by: 0)
                .map { id, meta1, cat1, meta2, cat2 -> [ id, meta1, cat1 + [cat2] ]}
                .set { all_cat }
        }

        // Process repeats and landscape analysis
        PROCESS_REPEATS(REPEATMASKER_WITH_OWN_LIB_2.out.masked, all_cat, Channel.fromPath( params.external_library, checkIfExists: true ))

        REPEATLANDSCAPE(PROCESS_REPEATS.out.align)

        DOWNLOAD_DFAM()

        REFORMAT_CLASSIFICATION_TO_DFAM(
            PROCESS_REPEATS.out.gff,
            PROCESS_REPEATS.out.cat,
            PROCESS_REPEATS.out.out,
            PROCESS_REPEATS.out.tbl,
            REPEATLANDSCAPE.out,
            DOWNLOAD_DFAM.out)

        REFORMAT_CLASSIFICATION_TO_DFAM.out.gff
            .join (genomes )
            .map { id, meta1, gff, meta2, fasta -> [ id, Utils.updateLinkedHashMap(meta2, 'repeats_gff',  meta2.repeats_gff + ['repeatmasker': gff]), fasta ] }
            .set { genomes }
        // GET_COMPLEX_REPEATS(CHANGE_CLASSIFICATION_TO_DFAM.out.gff)
        // SUMMARIZE_AND_PLOT_REPEATS(CHANGE_CLASSIFICATION_TO_DFAM.out.out.collect{it[1]}, newick)

    emit:
        genomes = genomes
    // masked_genomes = repeats_ch.fasta
}
