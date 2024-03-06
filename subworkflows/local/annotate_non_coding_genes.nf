workflow ANNOTATE_PROTEIN_CODING_GENES {
    take:
        genomes
        maker

    main:

    // select set about data to download
    RUN_INFERNAL(genomes)
    RUN_RNAMMER(genomes)
    RUN_BARNAP(genomes)

    // POST-PROCESSING: filter and concatenate results

    emit:
        annotated_genomes
}
