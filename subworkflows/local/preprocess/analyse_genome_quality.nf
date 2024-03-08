/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CALCULATE STATISTICS ABOUT GENOME QUALITY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// icnlude subworkflows
include { busco                             } from '../busco'

// include modules

workflow preprocess_genomes {
    take:
        genomes

    main:
    // check if there is key called busco_set in meta of genomes channel
    if ( ! genomes.meta[0].containsKey('busco_set')) {
        DOWNLOAD_BUSCO_DATASETS()

    }

    BUSCO(genomes_for_busco)
    REFORMAT_BUSCO(BUSCO.out.json)

    gather_genomes(REFORMAT_BUSCO.out)
        .set { genomes_for_busco_plotting }

    GATHER_BUSCO(genomes_for_busco_plotting.map{ ids, metas, files -> [ids, files, "busco", "csv", 'yes'] })
    GET_NEWICK_FOR_BUSCO(genomes_for_busco_plotting)

    // GATHER_BUSCO.out.view()
    // GET_NEWICK_FOR_BUSCO.out.newick.view()
    GET_NEWICK_FOR_BUSCO.out.newick
        .concat(GATHER_BUSCO.out)
        .map{it -> ["${it[0]}", it[1..-1]]}
        .groupTuple()
        .map{it -> [it[0], it[1][0][0], it[1][0][1], it[1][1]] }
        .set { genomes_for_busco_plotting }
    // GET_NEWICK_FOR_BUSCO.out.newick.view()
    // PLOT_BUSCO_SUMMARY(GATHER_BUSCO.out, )

    if (params.analysis_genome_quality ) {
        // add other statistics
    }

    emit:
    busco_on_assembly = REFORMAT_BUSCO.out
    all_busco_on_assembly = genomes_for_busco_plotting
}
