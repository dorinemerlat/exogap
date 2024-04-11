/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ANNOTATE PROTEIN CODING GENES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include subworkflows
include { PROCESS_MAKER as PROCESS_MAKER_BY_SIMILARITY_MAIN     } from './process_maker.nf'
include { PROCESS_MAKER as PROCESS_MAKER_BY_SIMILARITY_TRAINING } from './process_maker.nf'
include { PROCESS_MAKER as PROCESS_MAKER_AB_INITIO_1            } from './process_maker.nf'
include { PROCESS_MAKER as PROCESS_MAKER_AB_INITIO_2            } from './process_maker.nf'
include { CREATE_HMMS as CREATE_HMMS_1                          } from './create_hmms.nf'
include { CREATE_HMMS as CREATE_HMMS_2                          } from './create_hmms.nf'


// include modules
include { MAKER_BY_SIMILARITY as MAKER_BY_SIMILARITY_MAIN       } from '../../modules/local/maker_by_similarity'
include { MAKER_BY_SIMILARITY as MAKER_BY_SIMILARITY_TRAINING   } from '../../modules/local/maker_by_similarity'
include { MAKER_AB_INITIO as MAKER_AB_INITIO_1                  } from '../../modules/local/maker_ab_initio'
include { MAKER_AB_INITIO as MAKER_AB_INITIO_2                  } from '../../modules/local/maker_ab_initio'
include { FILTER_MAKER_AB_INITIO_PREDICTIONS                    } from '../../modules/local/filter_maker_ab_initio_predictions'
include { AGAT_MERGE_ANNOTATIONS                                } from '../../modules/local/agat_merge_annotations'
include { MAKER_MAP_IDS                                         } from '../../modules/local/maker_map_ids'
include { GENERATE_MAKER_FASTA                                  } from '../../modules/local/generate_maker_fasta'
include { BLASTP                                                } from '../../modules/local/blastp'
include { BLAST_FORMATTER as BLAST_FORMATTER_FMT2               } from '../../modules/local/blast_formatter'
include { BLAST_FORMATTER as BLAST_FORMATTER_FMT5               } from '../../modules/local/blast_formatter'
include { INTERPROSCAN                                          } from '../../modules/local/interproscan'
// include { FUNCTIONAL_ANNOTATION_BLAST2GO                            } from '../../modules/local/blast2go'
// include { FUNCTIONAL_ANNOTATION_WITH_MAKER                          } from '../../modules/local/functional_annotation_with_maker'


def gff_for_hmm( channel, iteration ) {
    return channel.map { id, meta, genome, gff, previous_iteration -> [ id, meta, genome, gff, iteration] }
}


def create
workflow ANNOTATE_PROTEIN_CODING_GENES {
    take:
        genomes
        blast_db

    main:
        // structural annotation
        genomes.map { id, meta, fasta -> [ id, meta, fasta, meta.main_protein_set.dataset.file, meta.transcript_set.dataset.file, meta.repeats_gff, 'main']}
            .set { genomes_for_main_similiarity_annotation }

        MAKER_BY_SIMILARITY_MAIN(genomes_for_main_similiarity_annotation)
        PROCESS_MAKER_BY_SIMILARITY_MAIN(MAKER_BY_SIMILARITY_MAIN.out)


        if (params.max_proteins_from_close_species ) {
            genomes.map { id, meta, fasta -> [ id, meta, fasta, meta.training_protein_set.dataset.file, meta.transcript_set.dataset.file, meta.repeats_gff, 'training']}
                .set { genomes_for_training_similiarity_annotation }

            MAKER_BY_SIMILARITY_TRAINING(genomes_for_training_similiarity_annotation)
            PROCESS_MAKER_BY_SIMILARITY_TRAINING(MAKER_BY_SIMILARITY_TRAINING.out)
            PROCESS_MAKER_BY_SIMILARITY_TRAINING.out.gff_for_augustus.set { gff_for_augustus_1 }

        }  else {
            PROCESS_MAKER_BY_SIMILARITY_MAIN.out.gff_for_augustus.set { gff_for_augustus_1 }
        }

        gff_for_snap_1 = gff_for_hmm(PROCESS_MAKER_BY_SIMILARITY_MAIN.out.gff_for_maker_and_snap, '1')
        gff_for_augustus_1 = gff_for_hmm(gff_for_augustus_1, '1')

        CREATE_HMMS_1(gff_for_snap_1, gff_for_augustus_1)

        PROCESS_MAKER_BY_SIMILARITY_MAIN.out.gff_for_maker_and_snap.join(CREATE_HMMS_1.out.snap)
            .map { id, meta1, genome1, gff, previous_iteration, meta2, genome2, snap, iteration -> [ id, meta1, genome1, gff, snap, iteration ] }
            .join(CREATE_HMMS_1.out.augustus)
            .map { id, meta1, genome1, gff, snap, iteration1, meta2, genome2, train, test, augustus, iteration2 -> [ id, meta1, genome1, gff, snap, augustus, iteration1 ] }
            .set { genomes_for_ab_initio_annotation_1 }

        MAKER_AB_INITIO_1(genomes_for_ab_initio_annotation_1)
        PROCESS_MAKER_AB_INITIO_1(MAKER_AB_INITIO_1.out)

        gff_for_snap_2 = gff_for_hmm(PROCESS_MAKER_AB_INITIO_1.out.gff_for_maker_and_snap, '2')
        gff_for_augustus_2 = gff_for_hmm(PROCESS_MAKER_AB_INITIO_1.out.gff_for_augustus, '2')

        FILTER_MAKER_AB_INITIO_PREDICTIONS(gff_for_augustus_2)
        FILTER_MAKER_AB_INITIO_PREDICTIONS.out.join(CREATE_HMMS_1.out.augustus_gff)
            .map { id, meta1, genome1, gff1, iteration1, meta2, genome2, gff2, iteration2 -> [ id, meta1, genome1, [gff1, gff2], iteration1 ] }
            .set { gff_for_augustus_2 }
        AGAT_MERGE_ANNOTATIONS(gff_for_augustus_2)

        CREATE_HMMS_2(gff_for_snap_2, AGAT_MERGE_ANNOTATIONS.out)

        PROCESS_MAKER_AB_INITIO_1.out.gff_for_maker_and_snap.join(CREATE_HMMS_2.out.snap)
            .map { id, meta1, genome1, gff, previous_iteration, meta2, genome2, snap, iteration -> [ id, meta1, genome1, gff, snap, iteration ] }
            .join(CREATE_HMMS_2.out.augustus)
            .map { id, meta1, genome1, gff, snap, iteration1, meta2, genome2, train, test, augustus, iteration2 -> [ id, meta1, genome1, gff, snap, augustus, iteration1 ] }
            .set { genomes_for_ab_initio_annotation_2 }

        MAKER_AB_INITIO_2(genomes_for_ab_initio_annotation_2)
        PROCESS_MAKER_AB_INITIO_2(MAKER_AB_INITIO_2.out)
        GENERATE_MAKER_FASTA(MAKER_AB_INITIO_2.out)

        PROCESS_MAKER_AB_INITIO_2.out.gff_for_maker_and_snap
            .join(GENERATE_MAKER_FASTA.out.proteins)
            .join(GENERATE_MAKER_FASTA.out.transcripts)
            .map { id, meta1, genome1, gff, iteration, meta2, genome2, proteins, meta3, genome3, transcripts -> [ id, meta1, genome1, gff, proteins, transcripts ]}
            .set { maker_outputs }

        MAKER_MAP_IDS(maker_outputs)

        // // functional annotation
        MAKER_MAP_IDS.out.proteins.combine(blast_db)
            .map { id, meta,  proteins, blastdb_name, blastdb_files -> [id, meta, proteins, blastdb_name, blastdb_files, '1.0E-6', '3', '20', 'functionnal']}
            .set { sequences_for_blastp }

        BLASTP(sequences_for_blastp)
        BLAST_FORMATTER_FMT2(BLASTP.out.map { id, meta, archive, comment -> [ id, meta, archive, '2', comment ]})
        BLAST_FORMATTER_FMT5(BLASTP.out.map { id, meta, archive, comment -> [ id, meta, archive, '5', comment ]})
        INTERPROSCAN(MAKER_MAP_IDS.out.proteins)

    //     if (params.blast2go) {
                // DOWNLOAD_OBO
                // DOWNLOAD_UNIPROT_BLAST()
    //         FUNCTIONAL_ANNOTATION_BLAST2GO(REFORMAT_MAKER_GFF, BLAST, INTERPROSCAN)
    //     } else {
    //         MAKER_FUNCTIONAL()
    //     }

    emit:
        gff = MAKER_MAP_IDS.out.gff
}
