/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ANNOTATE PROTEIN CODING GENES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include modules
include { SNAP                                      } from '../../modules/local/snap'
include { AGAT_SEPARATE_BY_RECORD_TYPE              } from '../../modules/local/agat_separate_by_record_type'
include { AGAT_FILTER_FEATURES_BY_ATTRIBUTE_VALUE   } from '../../modules/local/agat_filter_features_by_attribute_value'
include { AGAT_KEEP_LONGEST_ISOFORM                 } from '../../modules/local/agat_keep_longest_isoform'
include { AGAT_FILTER_BY_LOCUS_DISTANCE             } from '../../modules/local/agat_filter_by_locus_distance'
include { AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS } from '../../modules/local/agat_filter_incomplete_gene_coding_models'
include { AGAT_EXTRACT_SEQUENCE                     } from '../../modules/local/agat_extract_sequence'
include { MAKEBLASTDB                               } from '../../modules/local/makeblastdb'
include { BLASTP                                    } from '../../modules/local/blastp'
include { AGAT_FILTER_BY_MRNA_BLAST_VALUE           } from '../../modules/local/agat_filter_by_mrna_blast_value'
include { GFF_TO_AUGUSTUS_GFF                       } from '../../modules/local/gff_to_augustus_gff'
include { CREATE_AUGUSTUS_SET                       } from '../../modules/local/create_augustus_set'
include { CREATE_AUGUSTUS_NEW_SPECIES               } from '../../modules/local/create_augustus_new_species'
// include { TRAINING_AUGUSTUS                         } from '../../modules/local/training_augustus'
// include { OPTIMIZE_AUGUSTUS                         } from '../../modules/local/optimize_augustus'

workflow CREATE_HMMS {
    take:
        gff_for_snap
        gff_for_augustus

    main:
        // snap
        SNAP(gff_for_snap)

        // augustus
        // Separate according to type
        AGAT_SEPARATE_BY_RECORD_TYPE(gff_for_augustus)

        // Keep only those with a good aed
        AGAT_FILTER_FEATURES_BY_ATTRIBUTE_VALUE(AGAT_SEPARATE_BY_RECORD_TYPE.out.map { it << '0.1' })

        // Filter to keep longest isoforms
        AGAT_KEEP_LONGEST_ISOFORM(AGAT_FILTER_FEATURES_BY_ATTRIBUTE_VALUE.out)

        // Keep only genes with a certain distance between them
        AGAT_FILTER_BY_LOCUS_DISTANCE(AGAT_KEEP_LONGEST_ISOFORM.out)

        // Keep only complete genes
        AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS(AGAT_FILTER_BY_LOCUS_DISTANCE.out)

        // Eliminate redundancy
        AGAT_EXTRACT_SEQUENCE(AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS.out)
        MAKEBLASTDB(AGAT_EXTRACT_SEQUENCE.out.map {id, meta, genome, proteins, iteration -> [ id, meta, proteins, 'prot', iteration] })

        AGAT_EXTRACT_SEQUENCE.out.join(MAKEBLASTDB.out)
            .map { id, meta1, genome, proteins, iteration, meta2, db, comment -> [ id, meta1, proteins, db, '6', comment ] }
            .set { sequences_for_blastp }
        BLASTP(sequences_for_blastp)

        AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS.out.join(BLASTP.out)
            .map { id, meta1, genome, gff, iteration, meta2,  blastp_out, comment -> [ id, meta1, genome, gff, blastp_out, iteration ] }
            .set { agat_filter_by_mrna_blast_value_in }

        AGAT_FILTER_BY_MRNA_BLAST_VALUE(agat_filter_by_mrna_blast_value_in)

        // Reformat the gff file for augustus
        GFF_TO_AUGUSTUS_GFF(AGAT_FILTER_BY_MRNA_BLAST_VALUE.out)
        CREATE_AUGUSTUS_SET(GFF_TO_AUGUSTUS_GFF.out)

        CREATE_AUGUSTUS_NEW_SPECIES(gff_for_augustus)
        // CREATE_AUGUSTUS_SET.join(CREATE_AUGUSTUS_NEW_SPECIES.out)
        //     .map { -> [] }
        //     .set { training_set }

        // TRAINING_AUGUSTUS(training_set)
        // OPTIMIZE_AUGUSTUS(TRAINING_AUGUSTUS.out)
        // TRAINING_AUGUSTUS(OPTIMIZE_AUGUSTUS.out)

    emit:
        snap = SNAP.out
        // augustus = TRAINING_AUGUSTUS.out
}
