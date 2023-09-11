/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REFORMAT ASSEMBLY AND STATISTICS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTA_VALIDATOR } from '../../modules/local/fasta_validator'
include { GENOME_FORMAT } from '../../modules/local/genome_format'

workflow GENOME_PREPROCESS {
    take:
    genome_ch

    main:
    // Check if fasta file is valid
    FASTA_VALIDATOR(genome_ch)

    // Reformat genomes
    GENOME_FORMAT(FASTA_VALIDATOR.out)

    // Run some statistics on genome assemblies
    // GENOME_STATISTICS(GENOME_FORMAT.out.fasta)

    emit:
    fasta = GENOME_FORMAT.out.fasta
    // stats = GENOME_STATISTICS.out.stats
}
