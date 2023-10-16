/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REFORMAT ASSEMBLY AND STATISTICS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTA_VALIDATOR }                     from '../../modules/local/fasta_validator/run_fasta_validator'
include { RENAME_GENOMES }                      from '../../modules/local/bioawk/rename_genomes'
include { BUSCO }                               from '../../modules/local/busco/run_busco'
include { DOWNLOAD_BUSCO_DATASETS }             from '../../modules/local/busco/download_busco_datasets'
include { REFORMAT_BUSCO }                      from '../../modules/local/busco/reformat_busco'
include { GET_TAXONOMIC_LINEAGE }               from '../../modules/local/api/get_taxonomic-lineage'
include { GET_NEWICK as GET_NEWICK_FOR_ALL }    from '../../modules/local/api/get_newick'
include { GET_NEWICK as GET_NEWICK_FOR_BUSCO }  from '../../modules/local/api/get_newick'
include { CALCULATE_GENOME_SIZE }               from '../../modules/local/bioawk/calculate_genome_size'
include { GATHER_FILES as GATHER_BUSCO }        from '../../modules/local/gather-files'
include { PLOT_BUSCO_SUMMARY }                  from '../../modules/local/R_plots/plot_busco_summary'

def gather_genomes(genomes) {
    return genomes
        .map{it -> [it[0], it[1], it[2]]} // remove extra fiels
        .map { id, meta, file -> [id, ["${meta.taxid}": meta.name], file] } // keep only the IDs and essential meta data
        .flatten() // flat and list to get a destructed channel
        .toList()
        .map { genome -> genome.withIndex().collect { it, index -> [index % 3, it] } } // group by two: all the IDs, and all the meta data
        .flatMap { it } // have an item for IDs and one for meta data
        .groupTuple() // group to
        .map { index, it -> it.unique().sort() }
        .toList()
}

workflow PREPROCESS_GENOMES {
    take:
    genomes

    main:
    // Get the lineage taxonomy of each genome (one by genome)
    GET_TAXONOMIC_LINEAGE(genomes)

    GET_TAXONOMIC_LINEAGE.out
        .map{ it -> it.split('\n') }
        .map{ it -> [ it[0], it[1..-1].collect{ it.split(',') }]}
        .map{ id, lineage -> [ id, ["lineage": lineage.collect{ [ "rank": "${it[0]}", "name": "${it[2]}", "id": "${it[1]}"]}]]}
        .join(genomes)
        .map { id, taxonomy, meta, fasta -> [id, meta + taxonomy, fasta ] }
        .set { genomes }

    // Get the entire taxonomy tree (one for all genomes)
    GET_NEWICK_FOR_ALL(gather_genomes(genomes))

    // Check if fasta file is valid
    FASTA_VALIDATOR(genomes)

    // Reformat genomes
    RENAME_GENOMES(FASTA_VALIDATOR.out)

    RENAME_GENOMES.out.fasta
        .set { genomes }

    // Calculate the genome size
    CALCULATE_GENOME_SIZE(genomes)

    CALCULATE_GENOME_SIZE.out
        .map{id, meta, fasta, csv -> [ id, meta + ['assembly_size' : file(csv).readLines()[0].split(",")[1]], fasta ] }
        .set { genomes }

    // Busco on genome
    DOWNLOAD_BUSCO_DATASETS()

    genomes.map { id, meta, fasta -> meta.lineage.collect{[it.name.toLowerCase(), id, meta, fasta]} } // collect all parents of each genome
        .flatMap{ it }
        .combine(DOWNLOAD_BUSCO_DATASETS.out.splitCsv().map{ it -> [it[0].split('_')[0], it[0]] }, by:0) // combine with busco datasets
        .map{parent, id, meta, fasta, busco_set -> [id, meta, fasta, busco_set, 'genome']}
        .set{ genomes_for_busco }

    genomes_for_busco.groupTuple()
        .map{ it -> [ it[0], it[1][0] + ['busco_sets': it[3]], it[2][0] ]}
        .set{ genomes }

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

    // emit:
    // genomes     = genomes
    // newick      = GET_NEWICK.out.newick
}
