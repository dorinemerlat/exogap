/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CHECK IF FASTA FILES ARE VALIDS, REFORMATE THEM AND CALCULATE THEIR SIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNMASK_GENOME } from "$projectDir/modules/local/module_repeats/unmask_genome.nf"
include { REPEATMODELER } from "$projectDir/modules/local/module_repeats/repeatmodeler.nf"
include { HITE } from "$projectDir/modules/local/module_repeats/hite.nf"
include { CLEAN_REPEAT_FAMILIES } from "$projectDir/modules/local/module_repeats/clean_repeat_families.nf"
include { CD_HIT } from "$projectDir/modules/local/module_repeats/cd_hit.nf"
include { FAMILIES_CLUSTER_OVERLAP } from "$projectDir/modules/local/module_repeats/families_cluster_overlap.nf"
include { MCHELPER } from "$projectDir/modules/local/module_repeats/mchelper.nf"
include { GATHER_FAMILIES } from "$projectDir/modules/local/module_repeats/gather_families.nf"
include { CD_HIT as CD_HIT_FINAL } from "$projectDir/modules/local/module_repeats/cd_hit.nf"
include { REPEATMASKER } from "$projectDir/modules/local/module_repeats/repeatmasker.nf"

workflow REPEATS_ANNOTATION {
    take:
        genomes

    main:
        println "Run repeats annotation"

        // Some genomes are soft masked by default, it can be problematic for some tools
        UNMASK_GENOME(genomes)
        UNMASK_GENOME.out.set { genomes }

        // Run RepeatModeler to identify de novo repeats
        REPEATMODELER( genomes )

        // run HiTE
        HITE( genomes )

        // Clean up repeat libraries
        REPEATMODELER.out.map { id, meta, families -> [id, meta, families, file("$projectDir/data/repeats_classification.csv"), "RM2"] }
            .concat(HITE.out.confident_families.map { id, meta, families -> [id, meta, families, file("$projectDir/data/repeats_classification.csv"), "HC"] })
            .concat(HITE.out.low_confident_families.map { id, meta, families -> [id, meta, families, file("$projectDir/data/repeats_classification.csv"), "HLC"] })
            .set { repeat_families }

        CLEAN_REPEAT_FAMILIES( repeat_families )

        // // Merge libraries
        // CLEAN_REPEAT_FAMILIES.out.for_mchelper
        //     .groupTuple(by: [0,1])
        //     .map { id, meta, families_files -> [ id, meta, families_files, "nucleotide", "0.80", "10", "0" ] }
        //     // filter out where families_files contains not 4 files names
        //     .filter { id, meta, families_files, type, sequence_identity_threshold, length_of_throw_away_sequence, alignment_coverage_for_the_shorter_sequence ->
        //         families_files.size() == 3
        //     }
        //     .set { merged_libraries }

        // CD_HIT( merged_libraries )
        // FAMILIES_CLUSTER_OVERLAP(CD_HIT.out.clusters)

        // // Run MCHelper
        // merged_libraries.map { id, meta, families_files, type, sequence_identity_threshold, length_of_throw_away_sequence, alignment_coverage_for_the_shorter_sequence -> [id, meta, families_files] }
        //     .join(genomes)
        //     .map {id, meta1, families, meta2, genome -> [id, meta1, families, genome] }
        //     .set { all_families }

        // MCHELPER( all_families )

        // // Collect all results
        // MCHELPER.out.classified
        //     .join(MCHELPER.out.unclassified)
        //     .map { id, meta1, mchelper_classified, meta2, mchelper_unclassified -> [id, meta1, mchelper_classified, mchelper_unclassified]}
        //     .join(HITE.out.confident_families)
        //     .map { id, meta1, mchelper_classified, mchelper_unclassified, meta2, hite_confident_families -> [id, meta1, mchelper_classified, mchelper_unclassified, hite_confident_families]}
        //     .join(HITE.out.low_confident_families)
        //     .map { id, meta1, mchelper_classified, mchelper_unclassified, hite_confident_families, meta2, hite_low_confident_families -> [id, meta1, [mchelper_classified, mchelper_unclassified, hite_confident_families, hite_low_confident_families]]}
        //     .set { all_families_final }

        // GATHER_FAMILIES(all_families_final)

        // GATHER_FAMILIES.out
        //     .map { id, meta, families -> ['all', families] }
        //     .groupTuple(by: [0])
        //     .map{id, families_list -> [id, 'meta', [families_list, params.reference_library], "nucleotide", "0.80", "10", "0.8"]}
        //     .set { gathered_families }

        // CD_HIT_FINAL(gathered_families)

        // genomes
        //     .combine( CD_HIT_FINAL.out.library.map { id, meta, library -> library } )
        //     .set { genome_with_libraries }

        cd_hit_out = channel.fromPath( "/enadisk/tempor/merlat/exogaptwo/cache/repeats_annotation/cd_hit_final/all.nr.fa" ). map {file -> ['all', 'meta', file ] }

        genomes
            .combine( cd_hit_out.map { id, meta, library -> library } )
            .set { genome_with_libraries }

        // REPEATMASKER(genome_with_libraries)

        // GFF_TO_TSV(REPEATMASKER.out)
    // emit:
}

