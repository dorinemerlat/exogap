/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ANNOTATE PROTEIN CODING GENES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include subworkflows
include { AB_INITIO_ANNOTATION as AB_INITIO_ANNOTATION_1 }                  from 'ab_initio_annotation.nf'
include { AB_INITIO_ANNOTATION as AB_INITIO_ANNOTATION_2 }                  from 'ab_initio_annotation.nf'

// include modules
include { CTL_FOR_ANNOTATION_BY_SIMILARITY }                                from '../../modules/local/ctl_for_annotation_by_similarity'
include { MAKER_BY_SIMILARITY }                                             from '../../modules/local/maker_by_similarity'
include { GENERATE_MAKER_GFF }                                              from '../../modules/local/generate_maker_gff'
include { AGAT_STATISTICS }                                                 from '../../modules/local/agat_statistics'
include { ELIMINATE_REDUNDANCE as ELIMINATE_REDUNDANCE_IN_PROTEINS }        from '../../modules/local-hit/cd-hit_for_genes'
include { REFORMAT_MAKER_GFF }                                              from '../../modules/local/reformat_maker_gff'
include { BLAST }                                                           from '../../modules/local/blast'
include { INTERPROSCAN }                                                    from '../../modules/local/interproscan'
include { FUNCTIONAL_ANNOTATION_BLAST2GO }                                  from '../../modules/local/blast2go'
include { FUNCTIONAL_ANNOTATION_WITH_MAKER }                                from '../../modules/local/functional_annotation_with_maker'


workflow ANNOTATE_PROTEIN_CODING_GENES {
    take:
        genomes

    main:

        // structural annotation
        CTL_FOR_ANNOTATION_BY_SIMILARITY(genomes)
        MAKER_BY_SIMILARITY(CTL_FOR_ANNOTATION_BY_SIMILARITY.out)
        GENERATE_MAKER_GFF(MAKER_BY_SIMILARITY.out)
        AGAT_STATISTICS(GENERATE_MAKER_GFF.out)

        AB_INITIO_ANNOTATION_1(genomes, MAKER_BY_SIMILARITY)
        AB_INITIO_ANNOTATION_1(genomes, MAKER_AB_INITIO_1)
        REFORMAT_MAKER_GFF(MAKER_AB_INITIO_2.out)

        // functional annotation
        BLAST(REFORMAT_MAKER_GFF)
        INTERPROSCAN(REFORMAT_MAKER_GFF)

        if (params.blast2go) {
            FUNCTIONAL_ANNOTATION_BLAST2GO(REFORMAT_MAKER_GFF, BLAST, INTERPROSCAN)
        } else {
            FUNCTIONAL_ANNOTATION_WITH_MAKER(REFORMAT_MAKER_GFF, BLAST, INTERPROSCAN)
        }

    emit:
        annotated_genomes
}
