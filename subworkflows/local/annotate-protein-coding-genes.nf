workflow ANNOTATE_PROTEIN_CODING_GENES {
    take:
        masked_genomes

    main:
    MAKER_BY_SIMILARITY(masked_genomes)

    FILTER_MAKER_PREDICTIONS(MAKER_BY_SIMILARITY.out)
    HMM_WITH_AUGUSTUS_1(FILTER_MAKER_PREDICTIONS.out)
    HMM_WITH_SNAP_1(MAKER_BY_SIMILARITY.out)

    MAKER_AB_INITIO_1(HMM_WITH_AUGUSTUS_1.out, HMM_WITH_SNAP_1.out)

    FILTER_MAKER_PREDICTIONS(MAKER_AB_INITIO_1.out)
    HMM_WITH_AUGUSTUS_2(FILTER_MAKER_PREDICTIONS.out)
    HMM_WITH_SNAP_2(MAKER_BY_SIMILARITY.out)

    MAKER_AB_INITIO_2(HMM_WITH_AUGUSTUS_2.out, HMM_WITH_SNAP_2.out)

    REFORMAT_MAKER_GFF(MAKER_AB_INITIO_2.out)

    BLAST(REFORMAT_MAKER_GFF)
    INTERPROSCAN(REFORMAT_MAKER_GFF)
    BLAST2GO(REFORMAT_MAKER_GFF, BLAST, INTERPROSCAN)

    emit:
        annotated_genomes
}
