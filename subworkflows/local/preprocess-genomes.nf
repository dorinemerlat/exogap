/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REFORMAT ASSEMBLY AND STATISTICS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTA_VALIDATOR }                    from '../../modules/local/fasta-validator'
include { RENAME_GENOMES }                     from '../../modules/local/bioawk/rename-genomes'
include { BUSCO }                              from '../../modules/local/busco/busco'
include { DOWNLOAD_BUSCO_DATASETS }            from '../../modules/local/busco/download-busco-datasets'
include { SELECT_BUSCO_DATASETS }              from '../../modules/local/busco/select-busco-datasets'
include { GET_TAXONOMIC_LINEAGE }              from '../../modules/local/taxonomy/get-taxonomic-lineage'
include { GET_NEWICK }                         from '../../modules/local/taxonomy/get-newick'
include { CALCULATE_GENOME_SIZE }              from '../../modules/local/bioawk/calculate-genome-size'
include { GATHER_FILES as GATHER_GENOME_SIZE } from '../../modules/local/gather-files'

workflow PREPROCESS_GENOMES {
    take:
    genomes
    about_genomes

    main:
    // // Check if fasta file is valid
    // FASTA_VALIDATOR(genomes)

    // // Reformat genomes
    // RENAME_GENOMES(FASTA_VALIDATOR.out)

    // Calculate the genome size
    // CALCULATE_GENOME_SIZE(RENAME_GENOMES.out.fasta)
    // GATHER_GENOME_SIZE(CALCULATE_GENOME_SIZE.out.map{ it -> tuple('genome_sizes.tsv', it.last()) }.groupTuple())

    // // Get the genome taxonomy
    // GET_TAXONOMIC_LINEAGE(about_genomes)
    // taxonomy = GET_TAXONOMIC_LINEAGE.out.splitCsv(header: true)

    GET_NEWICK(about_genomes.map { it -> tuple( it.last()) }.collect(), about_genomes.toList())
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
