/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CHECK IF FASTA FILES ARE VALIDS, REFORMATE THEM AND CALCULATE THEIR SIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_LINEAGE }                    from '../../modules/local/download_lineage.nf'
include { DOWNLOAD_NEWICK }                     from '../../modules/local/download_newick'
include { FASTA_VALIDATOR }                     from '../../modules/local/fasta_validator'
include { RENAME_GENOME }                       from '../../modules/local/rename_genome'
include { CALCULATE_GENOME_SIZE }               from '../../modules/local/calculate_genome_size'

workflow PREPARE_GENOMES {
    take:
        genomes

    main:
        // Download the lineage and create newick file
        DOWNLOAD_LINEAGE(genomes)

        DOWNLOAD_LINEAGE.out
            .map { it -> [ it[0], file(it[1]).splitText() ] }
            .map{ id, lineage -> [ id, lineage.collect { it.split(',') } ] }
            .map{ id, lineage -> [ id, ["lineage": lineage.collect { [ "rank": "${it[0]}", "name": "${it[2]}".replace('\n', ''), "taxid": "${it[1]}" ] }]]}
            .join(genomes)
            .map { id, taxonomy, meta, fasta -> [id, meta + taxonomy, fasta ] }
            .set { genomes_with_lineage }

        DOWNLOAD_NEWICK(Utils.gatherGenomes(genomes_with_lineage))

        // Check if fasta file is valid
        // FASTA_VALIDATOR(genomes)

        // Reformat genomes
        RENAME_GENOME(genomes_with_lineage)

        // Calculate the genome size
        CALCULATE_GENOME_SIZE(RENAME_GENOME.out.fasta)
        CALCULATE_GENOME_SIZE.out
            .map{id, meta, fasta, csv -> [ id, meta + ['assembly_size' : file(csv).readLines()[0].split(",")[1]], fasta ] }
            .map{id, meta, fasta -> [ id, meta + ['genome' : fasta], fasta ] }
            .set { genomes_with_info }

    emit:
        genomes     = genomes_with_info
}
