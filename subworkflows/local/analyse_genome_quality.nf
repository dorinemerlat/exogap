/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CALCULATE STATISTICS ABOUT GENOME QUALITY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// include subworkflows
include { BUSCO                             } from '../../modules/local/busco'
include { DOWNLOAD_BUSCO_DATASETS           } from '../../modules/local/download_busco_datasets'
include { GATHER_FILES as GATHER_BUSCO      } from '../../modules/local/gather_files'
// include modules

workflow ANALYSE_GENOME_QUALITY {
    take:
        genomes 
        newick
        info
        
    main:
        // choose which busco dataset to use for each specie
        DOWNLOAD_BUSCO_DATASETS()

        genomes.map { id, meta, fasta -> meta.lineage.collect{[it.name.toLowerCase(), id, meta, fasta]} } // create an it for each parent of each genome
            .flatMap{ it }
            .combine(DOWNLOAD_BUSCO_DATASETS.out.splitCsv(), by: 0) // -> [taxon, id, meta, fasta, busco_dataset]
            .map { taxon, id, meta, fasta, busco_dataset -> [id, meta, fasta, busco_dataset] }
            .set { genomes_for_busco }

        genomes_for_busco.groupTuple()
            .map { id, meta, file, busco_datasets -> [id, meta[0] + [ "busco_dataset": busco_datasets ], file[0]] }
            .set { genomes }

        BUSCO(genomes_for_busco.map {id, meta, fasta, busco_dataset -> [id, meta, fasta, busco_dataset, 'genome' ]} )
        BUSCO.out.csv
            .groupTuple()
            .map { id, files -> [id, "meta", files, "busco.csv", "yes"] }
            .set { busco_to_gather }

        GATHER_FILES(busco_to_gather)
        JOIN_FILES(info.combine( GATHER_BUSCO.out.map {dataset, meta, busco -> [dataset, busco] }).map { info, dataset, busco -> ['info_and_busco_' + dataset, info, busco, "info.csv"] })
        // PLOT_BUSCO(GATHER_FILES.out, newick, info)
    //     REFORMAT_BUSCO(BUSCO.out.json)

    // gather_genomes(REFORMAT_BUSCO.out)
    //     .set { genomes_for_busco_plotting }

    // GATHER_BUSCO(genomes_for_busco_plotting.map{ ids, metas, files -> [ids, files, "busco", "csv", 'yes'] })
    // GET_NEWICK_FOR_BUSCO(genomes_for_busco_plotting)

    // // GATHER_BUSCO.out.view()
    // // GET_NEWICK_FOR_BUSCO.out.newick.view()
    // GET_NEWICK_FOR_BUSCO.out.newick
    //     .concat(GATHER_BUSCO.out)
    //     .map{it -> ["${it[0]}", it[1..-1]]}
    //     .groupTuple()
    //     .map{it -> [it[0], it[1][0][0], it[1][0][1], it[1][1]] }
    //     .set { genomes_for_busco_plotting }
    // // GET_NEWICK_FOR_BUSCO.out.newick.view()
    // // PLOT_BUSCO_SUMMARY(GATHER_BUSCO.out, )

    // if (params.analysis_genome_quality ) {
    //     // add other statistics
    // }

    // emit:
    // busco_on_assembly = REFORMAT_BUSCO.out
    // all_busco_on_assembly = genomes_for_busco_plotting
}
