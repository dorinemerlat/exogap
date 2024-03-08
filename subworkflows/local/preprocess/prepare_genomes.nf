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
        fasta_validator(genomes)

        // Reformat genomes
        rename_genome(fasta_validator.out)

        rename_genome.out.fasta
            .set { genomes }

        // Calculate the genome size
        calculate_genome_size(genomes)

        calculate_genome_size.out
            .map{id, meta, fasta, csv -> [ id, meta + ['assembly_size' : file(csv).readLines()[0].split(",")[1]], fasta ] }
            .set { genomes }


    emit:
        genomes     = genomes
}
