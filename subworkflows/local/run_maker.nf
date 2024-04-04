/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RUN MAKER AND PROCESS RESULTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// include modules
include { MAKER                            } from '../../modules/local/maker'
include { GENERATE_MAKER_GFF               } from '../../modules/local/generate_maker_gff'
include { AGAT_STATISTICS                  } from '../../modules/local/agat_statistics'
include { GENERATE_AED                     } from '../../modules/local/generate_aed'

workflow RUN_MAKER {
    take:
        genomes

    main:
        MAKER(genomes)
        GENERATE_MAKER_GFF(MAKER.out)
        AGAT_STATISTICS(GENERATE_MAKER_GFF.out.gff)
        GENERATE_AED(GENERATE_MAKER_GFF.out.gff)

    emit:
        gff_for_snap = GENERATE_MAKER_GFF.out.gff_complete_with_sequences
        gff_for_augustus = GENERATE_MAKER_GFF.out.gff
    //     stats = AGAT_STATISTICS.out
}
