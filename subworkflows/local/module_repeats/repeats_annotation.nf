/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CHECK IF FASTA FILES ARE VALIDS, REFORMATE THEM AND CALCULATE THEIR SIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNMASK_GENOME } from "$projectDir/modules/local/module_repeats/unmask_genome.nf"
include { REPEATMODELER } from "$projectDir/modules/local/module_repeats/repeatmodeler.nf"
include { HITE } from "$projectDir/modules/local/module_repeats/hite.nf"
// include { EARLGREY } from "$projectDir/modules/local/module_repeats/earlgrey.nf"
// include { MCHELPER } from "$projectDir/modules/local/module_repeats/mchelper.nf"

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

        // get dfam files: if famdb_path is a directory, collect all .h5 files inside it
        // if (params.famdb_partitions)  {
        //     println "Using Dfam database files: ${params.famdb_partitions}"
        // }
        // // Run EarlGrey to filter and better classify the repeats
        // genomes
        //     .join( REPEATMODELER.out )
        //     .map { id, meta1, genome, meta2, families -> [id, meta1, genome, families, params.famdb_partitions] }
        //     .set { genomes_and_families }

        // EARLGREY( genomes_and_families )

        // run HiTE
        HITE( genomes )

        // // Run MCHelper
        // MCHELPER( EARLGREY.out, genomes )

    // emit:
}

