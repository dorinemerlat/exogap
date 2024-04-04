/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ANNOTATE PROTEIN CODING GENES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include subworkflows
include { RUN_MAKER as RUN_MAKER_BY_ANNOTATION_MAIN     } from './run_maker.nf'
include { RUN_MAKER as RUN_MAKER_BY_ANNOTATION_TRAINING } from './run_maker.nf'
include { RUN_MAKER as RUN_AB_INITIO_MAKER_1            } from './run_maker.nf'
include { RUN_MAKER as RUN_AB_INITIO_MAKER_2            } from './run_maker.nf'
include { CREATE_HMMS as CREATE_HMMS_1                  } from './create_hmms.nf'
include { CREATE_HMMS as CREATE_HMMS_2                  } from './create_hmms.nf'


// include modules
include { CTL_FOR_ANNOTATION_BY_SIMILARITY as CTL_FOR_ANNOTATION_BY_SIMILARITY_MAIN     } from '../../modules/local/ctl_for_annotation_by_similarity'
include { CTL_FOR_ANNOTATION_BY_SIMILARITY as CTL_FOR_ANNOTATION_BY_SIMILARITY_TRAINING } from '../../modules/local/ctl_for_annotation_by_similarity'
include { CTL_FOR_AB_INITIO_ANNOTATION as CTL_FOR_AB_INITIO_ANNOTATION_1                } from '../../modules/local/ctl_for_ab_initio_annotation'
include { CTL_FOR_AB_INITIO_ANNOTATION as CTL_FOR_AB_INITIO_ANNOTATION_2                } from '../../modules/local/ctl_for_ab_initio_annotation'

include { REFORMAT_MAKER_GFF                                        } from '../../modules/local/reformat_maker_gff'
include { BLASTP                                                    } from '../../modules/local/blastp'
// include { INTERPROSCAN                                              } from '../../modules/local/interproscan'
// include { FUNCTIONAL_ANNOTATION_BLAST2GO                            } from '../../modules/local/blast2go'
// include { FUNCTIONAL_ANNOTATION_WITH_MAKER                          } from '../../modules/local/functional_annotation_with_maker'


workflow ANNOTATE_PROTEIN_CODING_GENES {
    take:
        genomes

    main:

        // structural annotation
        genomes.map { id, meta, fasta -> [ id, meta, fasta, meta.main_protein_set.dataset.file, meta.transcript_set.dataset.file, meta.repeats_gff, 'main']}
            .set { genomes_for_main_similiarity_annotation }

        CTL_FOR_ANNOTATION_BY_SIMILARITY_MAIN(genomes_for_main_similiarity_annotation)
        RUN_MAKER_BY_ANNOTATION_MAIN(CTL_FOR_ANNOTATION_BY_SIMILARITY_MAIN.out)


        if (params.max_proteins_from_close_species ) {
            genomes.map { id, meta, fasta -> [ id, meta, fasta, meta.training_protein_set.dataset.file, meta.transcript_set.dataset.file, meta.repeats_gff, 'training']}
                .set { genomes_for_training_similiarity_annotation }

            CTL_FOR_ANNOTATION_BY_SIMILARITY_TRAINING(genomes_for_training_similiarity_annotation)
            RUN_MAKER_BY_ANNOTATION_TRAINING(CTL_FOR_ANNOTATION_BY_SIMILARITY_TRAINING.out)
            RUN_MAKER_BY_ANNOTATION_TRAINING.out.gff_for_augustus.set { gff_for_augustus }

        }  else {
            RUN_MAKER_BY_ANNOTATION_MAIN.out.gff_for_augustus.set { gff_for_augustus }
        }

        RUN_MAKER_BY_ANNOTATION_MAIN.out.gff_for_snap
            .map { id, meta, genome, gff, previous_iteration -> [ id, meta, genome, gff, '1'] }
            .set { gff_for_snap }

        gff_for_augustus
            .map { id, meta, genome, gff, previous_iteration -> [ id, meta, genome, gff, '1'] }
            .set { gff_for_augustus }

        CREATE_HMMS_1(gff_for_snap, gff_for_augustus)
    //     CTL_FOR_AB_INITIO_1(RUN_MAKER_BY_ANNOTATION_MAIN.out, CREATE_HMMS_1.out.snap, CREATE_HMMS_1.out.augustus)
    //     RUN_AB_INITIO_MAKER_1(CTL_FOR_AB_INITIO_1.out, '1')

    //     // TO DO: create set for augustus2
    //     CREATE_HMMS_2(RUN_AB_INITIO_MAKER_1.out, agustus_set2, '2')
    //     CTL_FOR_AB_INITIO_2(RUN_AB_INITIO_MAKER_1.out, CREATE_HMMS_2.out.snap, CREATE_HMMS_2.out.augustus, '2')
    //     RUN_AB_INITIO_MAKER_2(CTL_FOR_AB_INITIO_2.out, '2')

    //     // functional annotation
    //     REFORMAT_PROTEINS(RUN_AB_INITIO_MAKER_2.out) // TO DO: add this to the module to remove star
    //     BLASTP(REFORMAT_MAKER_GFF)
    //     INTERPROSCAN(REFORMAT_MAKER_GFF)

    //     if (params.blast2go) {
                // DOWNLOAD_OBO
                // DOWNLOAD_UNIPROT_BLAST()
    //         FUNCTIONAL_ANNOTATION_BLAST2GO(REFORMAT_MAKER_GFF, BLAST, INTERPROSCAN)
    //     } else {
    //         FUNCTIONAL_ANNOTATION_WITH_MAKER(REFORMAT_MAKER_GFF, BLAST, INTERPROSCAN)
    //     }

    // emit:
    //     annotated_genomes
}
