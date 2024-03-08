/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ANNOTATE PROTEIN CODING GENES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include subworkflows
include { AB_INITIO_ANNOTATION as AB_INITIO_ANNOTATION_1 }                  from 'ab_initio_annotation.nf'
include { AB_INITIO_ANNOTATION as AB_INITIO_ANNOTATION_2 }                  from 'ab_initio_annotation.nf'

// include modules
include { MAKER_BY_SIMILARITY }                                             from '../../modules/local/maker/maker_by_similarity'
include { ELIMINATE_REDUNDANCE as ELIMINATE_REDUNDANCE_IN_PROTEINS }        from '../../modules/local/cd-hit/cd-hit_for_genes'
include { REFORMAT_MAKER_GFF }                                              from '../../modules/local/maker/reformat_maker_gff'
include { BLAST }                                                           from '../../modules/local/blast/blast'
include { INTERPROSCAN }                                                    from '../../modules/local/interproscan/interproscan'
include { FUNCTIONAL_ANNOTATION_BLAST2GO }                                  from '../../modules/local/blast2go/blast2go'
include { FUNCTIONAL_ANNOTATION_WITH_MAKER }                                from '../../modules/local/maker/functional_annotation_with_maker'


workflow ANNOTATE_PROTEIN_CODING_GENES {
    take:
        genomes

    main:

        // download protein public data
        DOWNLOAD_PROTEINS_FROM_PROTEINS(SELECT_PROTEINS_FROM_PROTEOMES.out)
        DOWNLOAD_PROTEINS_FROM_CLOSEST_PROTEINS(SELECT_PROTEINS_FROM_CLOSEST_SPECIES.out)

        ELIMINATE_REDUNDANCE_IN_PROTEINS(DOWNLOAD_PROTEINS.out, DOWNLOAD_PROTEINS.out)

        // download transcript public data
        DOWNLOAD_TSA_TRANSCRIPTS(SELECT_TRANSCRIPTS_FROM_TSA.out)
        DOWNLOAD_SRA_TRANSCRIPTS(SELECT_TRANSCRIPTS_FROM_SRA.out)
        ASSEMBLE_SRA_READS(DOWNLOAD_SRA_TRANSCRIPTS.out)
        ELIMINATE_REDUNDANCE_IN_TRANSCRIPTS(DOWNLOAD_TSA_TRANSCRIPTS.out, DOWNLOAD_SRA_TRANSCRIPTS.out)

        // structural annotation
        MAKER_BY_SIMILARITY(genomes, ELIMINATE_REDUNDANCE_IN_PROTEINS, ELIMINATE_REDUNDANCE_IN_TRANSCRIPTS)
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
