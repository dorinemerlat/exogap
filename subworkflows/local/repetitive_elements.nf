include {REPEATMODELER          } from '../../modules/local/repeatmodeler'

workflow REPETITIVE_ELEMENTS{
    take:
        ch_genomes

    main:
    REPEATMODELER(ch_genomes)

    emit:
    repeats = REPEATMODELER.out
}

