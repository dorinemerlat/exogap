/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CHECK IF FASTA FILES ARE VALIDS, REFORMATE THEM AND CALCULATE THEIR SIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { fasta_validator }                     from '../../../modules/local/fasta_validator/fasta_validator'
include { rename_genome }                       from '../../../modules/local/bioawk/rename_genome'
include { calculate_genome_size }               from '../../../modules/local/bioawk/calculate_genome_size'

workflow prepare_genomes {
    take:
        genomes

    main:

        // Check if fasta file is valid
        FASTA_VALIDATOR(genomes)

        // Reformat genomes
        RENAME_GENOME(FASTA_VALIDATOR.out)

        RENAME_GENOME.out.fasta
            .set { genomes }

        // Calculate the genome size
        CALCULATE_GENOME_SIZE(genomes)

        CALCULATE_GENOME_SIZE.out
            .map{id, meta, fasta, csv -> [ id, meta + ['assembly_size' : file(csv).readLines()[0].split(",")[1]], fasta ] }
            .set { genomes }


    emit:
        genomes     = genomes
}
