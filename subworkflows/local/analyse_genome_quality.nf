/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CALCULATE STATISTICS ABOUT GENOME QUALITY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// include subworkflows
include { BUSCO                             } from '../../modules/local/busco'
include { DOWNLOAD_BUSCO_LIST_DATASETS      } from '../../modules/local/download_busco_list_datasets'
include { DOWNLOAD_BUSCO_DATASETS           } from '../../modules/local/download_busco_datasets'
include { GATHER_FILES as GATHER_BUSCO      } from '../../modules/local/gather_files'
include { BUSCO_TO_CSV                      } from '../../modules/local/busco_to_csv'  
include { JOIN_FILES                        } from '../../modules/local/join_files'
// include { PLOT_BUSCO                        } from '../../modules/local/plot_busco'

workflow ANALYSE_GENOME_QUALITY {
    take:
        genomes 
        newick
        info
        
    main:
        // choose which busco dataset to use for each specie
        DOWNLOAD_BUSCO_LIST_DATASETS()
        genomes.map { id, meta, fasta -> meta.lineage.collect{[it.name.toLowerCase(), id, meta, fasta]} } // create an it for each parent of each genome
            .flatMap{ it }
            .combine(DOWNLOAD_BUSCO_LIST_DATASETS.out.splitCsv(), by: 0) // -> [taxon, id, meta, fasta, busco_dataset]
            .map { taxon, id, meta, fasta, busco_dataset -> [id, meta, fasta, busco_dataset] }
            .set { genomes_for_busco }

        genomes_for_busco.groupTuple()
            .map { id, meta, fasta, busco_datasets -> [id, meta[0] + [ "busco_dataset": busco_datasets ], fasta[0]] }
            .set { genomes }
        
        genomes_for_busco. map { id, meta, fasta, busco_dataset -> busco_dataset }
            .unique()
            .collect()
            // join the list it to a string with space as separator
            .map {it -> it.join(" ")}
            .map { it -> ['busco_dataset', it]}
            .set { busco_datasets }

        DOWNLOAD_BUSCO_DATASETS(busco_datasets)
        
        genomes_for_busco.combine(DOWNLOAD_BUSCO_DATASETS.out)
            .map {id, meta, fasta, dataset, busco_downloads -> [id, meta, fasta, dataset, 'genome', busco_downloads] }
            .set { genomes_for_busco }
        
        BUSCO(genomes_for_busco)
        BUSCO_TO_CSV(BUSCO.out.json)
        BUSCO_TO_CSV.out
            .groupTuple()
            .map { id, files -> [id, "meta", files, "busco.csv", "yes"] }
            .set { busco_to_gather }

        GATHER_BUSCO(busco_to_gather)
        info.combine( GATHER_BUSCO.out.map {dataset, meta, busco -> [dataset, busco] }).map { info, dataset, busco -> ['info_and_busco_' + dataset, [info, busco], "info.csv"] }.view()
        // JOIN_FILES(info.combine( GATHER_BUSCO.out.map {dataset, meta, busco -> [dataset, busco] }).map { info, dataset, busco -> ['info_and_busco_' + dataset, [info, busco], "info.csv"] })
        // PLOT_BUSCO(GATHER_FILES.out, newick, info)

    // emit:
    //     busco_datasets = DOWNLOAD_BUSCO_DATASETS.out
    // busco_on_assembly = REFORMAT_BUSCO.out
    // all_busco_on_assembly = genomes_for_busco_plotting
}
