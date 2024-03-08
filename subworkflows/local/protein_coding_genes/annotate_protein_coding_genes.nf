/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ANNOTATE PROTEIN CODING GENES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include subworkflows
include { ab_initio_annotation as ab_initio_annotation_1 }                  from 'ab_initio_annotation.nf'
include { ab_initio_annotation as ab_initio_annotation_2 }                  from 'ab_initio_annotation.nf'

// include modules
include { download_proteins as download_proteins_from_proteins }            from '../../modules/local/api/download_proteins'
include { download_proteins as download_proteins_from_closest_proteins }    from '../../modules/local/api/download_proteins'
include { download_tsa_transcripts }                                        from '../../modules/local/api/download_tsa_transcripts'
include { download_sra_transcripts }                                        from '../../modules/local/api/download_sra_transcripts'
include { assemble_sra_reads }                                              from '../../modules/local/star/assemble_sra_reads'
include { eliminate_redundance as eliminate_redundance_in_proteins }        from '../../modules/local/cd-hit/cd-hit_for_genes'
include { eliminate_redundance as eliminate_redundance_in_transcripts }     from '../../modules/local/cd-hit/cd-hit-est'
include { maker_by_similarity }                                             from '../../modules/local/maker/maker_by_similarity'
include { eliminate_redundance as eliminate_redundance_in_proteins }        from '../../modules/local/cd-hit/cd-hit_for_genes'
include { reformat_maker_gff }                                              from '../../modules/local/maker/reformat_maker_gff'
include { blast }                                                           from '../../modules/local/blast/blast'
include { interproscan }                                                    from '../../modules/local/interproscan/interproscan'
include { functional_annotation_blast2go }                                  from '../../modules/local/blast2go/blast2go'
include { functional_annotation_with_maker }                                from '../../modules/local/maker/functional_annotation_with_maker'


workflow annotate_protein_coding_genes {
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
