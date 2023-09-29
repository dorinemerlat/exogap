/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REFORMAT ASSEMBLY AND STATISTICS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTA_VALIDATOR }                     from '../../modules/local/fasta-validator'
include { RENAME_GENOMES }                      from '../../modules/local/bioawk/rename-genomes'
include { BUSCO }                               from '../../modules/local/busco/busco'
include { DOWNLOAD_BUSCO_DATASETS }             from '../../modules/local/busco/download-busco-datasets'
include { SELECT_BUSCO_DATASETS }               from '../../modules/local/busco/select-busco-datasets'
include { GET_TAXONOMIC_LINEAGE }               from '../../modules/local/taxonomy/get-taxonomic-lineage'
include { GET_NEWICK }                          from '../../modules/local/taxonomy/get-newick'
include { CALCULATE_GENOME_SIZE }               from '../../modules/local/bioawk/calculate-genome-size'
include { GATHER_FILES as GATHER_GENOME_SIZE }  from '../../modules/local/gather-files'

workflow PREPROCESS_GENOMES {
    take:
    genomes

    main:
    // Get the genome taxonomy
    GET_TAXONOMIC_LINEAGE(genomes)

    genomes = GET_TAXONOMIC_LINEAGE.out
                .splitJson()
                .buffer( size: 11 )
                .map{ it -> tuple( file("${it[0].value}"), [ 'lineage': tuple(
                        ["${it[1].key}" :it[1].value ],
                        ["${it[2].key}" :it[2].value ],
                        ["${it[3].key}" :it[3].value ],
                        ["${it[4].key}" :it[4].value ],
                        ["${it[5].key}" :it[5].value ],
                        ["${it[6].key}" :it[6].value ],
                        ["${it[7].key}" :it[7].value ],
                        ["${it[8].key}" :it[8].value ],
                        ["${it[9].key}" :it[9].value ],
                        ["${it[10].key}" :it[1].value ],
                        )])}
                .join(genomes.map { meta, it -> tuple(file(it), meta) })
                .map { fasta, taxonomy, meta -> [meta + taxonomy , file(fasta) ] }
                .view()

    // GET_NEWICK(genomes.map{ meta, file -> [ "${meta.taxid}": meta.name ] }.collect(flat : false))

    // // Check if fasta file is valid
    // FASTA_VALIDATOR(genomes)

    // // Reformat genomes
    // RENAME_GENOMES(FASTA_VALIDATOR.out)

    // // Calculate the genome size
    // CALCULATE_GENOME_SIZE(RENAME_GENOMES.out.fasta)
    // GATHER_GENOME_SIZE(CALCULATE_GENOME_SIZE.out.map{ it -> tuple('genome_sizes.tsv', it.last()) }.groupTuple())

    // // Get the genome taxonomy
    // GET_TAXONOMIC_LINEAGE(about_genomes)
    // taxonomy = GET_TAXONOMIC_LINEAGE.out.splitCsv(header: true)

    // GET_NEWICK(about_genomes.map { it -> tuple( it.last()) }.collect(), about_genomes.toList())
    // // // Busco on genome
    // // DOWNLOAD_BUSCO_DATASETS()
    // // SELECT_BUSCO_DATASETS(DOWNLOAD_BUSCO_DATASETS.out, GET_TAXONOMIC_LINEAGE.out)
    // // BUSCO(RENAME_GENOMES.out.fasta, SELECT_BUSCO_DATASETS.out.splitCsv(), 'genome')


    // emit:
    // fasta = RENAME_GENOMES.out.fasta
    // size = CALCULATE_GENOME_SIZE.out
    // taxonomy = taxonomy
    // taxonomy = DOWNLOAD_BUSCO_DATASETS.out
    // busco = SELECT_BUSCO_DATASETS
    // stats = GENOME_STATISTICS.out.stats
}
