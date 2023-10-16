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
include {REFORMAT_OUT                                               } from '../../modules/local/repeatmasker/reformat-out'
include {SUMMARIZE_REPEATS                                          } from '../../modules/local/stats-and-plots/summarize-repeats'
// include {STATS                                                      } from '../../modules/local/re_stats/stats'
// include {PLOTS                                                      } from '../../modules/local/re_stats/plots'

workflow REPETITIVE_ELEMENTS {
    take:
        genomes
        newick

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

        // // make 2 iterations with own library (all or specie-specific)
        // REPEATMASKER_WITH_OWN_LIB_1(REPEATMASKER_WITH_EXISTING_LIB.out.masked, CD_HIT_FOR_REPEATS_1.out)
        // REPEATMASKER_WITH_OWN_LIB_2(REPEATMASKER_WITH_OWN_LIB_1.out.masked, CD_HIT_FOR_REPEATS_2.out)

        // // Group all .cat files
        // REPEATMASKER_WITH_EXISTING_LIB.out.cat
        //     .concat(REPEATMASKER_WITH_OWN_LIB_1.out.cat, REPEATMASKER_WITH_OWN_LIB_2.out.cat)
        //     .groupTuple()
        //     .set { cat }

        // // Process repeats and landscape analysis
        // PROCESS_REPEATS(REPEATMASKER_WITH_OWN_LIB_2.masked, cat, external_library)
        // REPEATLANDSCAPE(PROCESS_REPEATS.out.align)

/////////
        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .set {genomes_for_join}

        // Masked
        Channel.fromPath('jolly_watson.log.workdir').map {it -> file(it).readLines()}.flatten().map {it -> [file("${it}*.align").simpleName, it]}
            .map {name, it -> tuple(name[0], file("${it}${name[0]}.masked")) }
            .set {masked}
        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .join(masked).map { id, meta, genome, file -> [ meta, file ] }
            .set {masked}

        // gff
        Channel.fromPath('jolly_watson.log.workdir').map {it -> file(it).readLines()}.flatten().map {it -> [file("${it}*.align").simpleName, it]}
            .map {name, it -> tuple(name[0], file("${it}${name[0]}.gff")) }
            .set {gff}
        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .join(gff).map { id, meta, genome, file -> [ meta, file ] }
            .set {gff}

        // complex gff
        Channel.fromPath('jolly_watson.log.workdir').map {it -> file(it).readLines()}.flatten().map {it -> [file("${it}*.align").simpleName, it]}
            .map {name, it -> tuple(name[0], file("${it}${name[0]}.-complex.gff")) }
            .set {complex_gff}
        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .join(complex_gff).map { id, meta, genome, file -> [ meta, file ] }
            .set {complex_gff}

        // align
        Channel.fromPath('jolly_watson.log.workdir').map {it -> file(it).readLines()}.flatten().map {it -> [file("${it}*.align").simpleName, it]}
            .map {name, it -> tuple(name[0], file("${it}${name[0]}.align")) }
            .set {align}
        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .join(align).map { id, meta, genome, file -> [ meta, file ] }
            .set {align}

        // cat
        Channel.fromPath('jolly_watson.log.workdir').map {it -> file(it).readLines()}.flatten().map {it -> [file("${it}*.align").simpleName, it]}
            .map {name, it -> tuple(name[0], file("${it}${name[0]}.cat")) }
            .set {cat}
        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .join(cat).map { id, meta, genome, file -> [ meta, file ] }
            .set {cat}

        //out
        Channel.fromPath('jolly_watson.log.workdir').map {it -> file(it).readLines()}.flatten().map {it -> [file("${it}*.align").simpleName, it]}
            .map {name, it -> tuple(name[0], file("${it}${name[0]}.out")) }
            .set {out}
        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .join(out).map { id, meta, genome, file -> [ meta, file ] }
            .set {out}

        //tbl
        Channel.fromPath('jolly_watson.log.workdir').map {it -> file(it).readLines()}.flatten().map {it -> [file("${it}*.align").simpleName, it]}
            .map {name, it -> tuple(name[0], file("${it}${name[0]}.tbl")) }
            .set {tbl}
        genomes.map { meta, file -> [ meta.id, meta, file ] }
            .join(tbl).map { id, meta, genome, file -> [ meta, file ] }
            .set {tbl}

//////////
        Channel.fromPath('data/repetitive_elements/TEClasses.tsv').set {classification}
        REFORMAT_OUT(out)

        REFORMAT_OUT.out.out.map {meta, it -> [it]}.collect()

        DOWNLOAD_TE_CLASSES()

        SUMMARIZE_REPEATS(REFORMAT_OUT.out.out.map {meta, it -> [it]}.collect(), classification, newick)

        // REPEATLANDSCAPE(align)

        // EXTRACT_SEQUENCES(masked, gff)

        // PLOT_REPEATLANDSCAPE(REPEATLANDSCAPE.out) -> one by genome
        // PlOT_REPEAT_OVERVIEW() -> one for all
        // PLOT_CORRELATION() -> one for each group
        // PLOT_ABUNDANCE() -> one for all

    // emit:
    // re_gff3 = repeats_ch.gff3
    // masked_genomes = repeats_ch.fasta
}
