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
        download_busco_datasets()

    }

    BUSCO(genomes_for_busco)
    reformat_busco(busco.out.json)

    gather_genomes(reformat_busco.out)
        .set { genomes_for_busco_plotting }

    GATHER_BUSCO(genomes_for_busco_plotting.map{ ids, metas, files -> [ids, files, "busco", "csv", 'yes'] })
    GET_NEWICK_FOR_BUSCO(genomes_for_busco_plotting)

    // gather_busco.out.view()
    // get_newick_for_busco.out.newick.view()
    get_newick_for_busco.out.newick
        .concat(gather_busco.out)
        .map{it -> ["${it[0]}", it[1..-1]]}
        .groupTuple()
        .map{it -> [it[0], it[1][0][0], it[1][0][1], it[1][1]] }
        .set { genomes_for_busco_plotting }
    // get_newick_for_busco.out.newick.view()
    // plot_busco_summary(gather_busco.out, )

    if (params.analysis_genome_quality ) {
        // add other statistics
    }

    emit:
    busco_on_assembly = reformat_busco.out
    all_busco_on_assembly = genomes_for_busco_plotting
}
