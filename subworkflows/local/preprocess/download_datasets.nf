/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DOWNLOAD SELECTED PUBLIC DATA TO USE FOR ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include modules

// include {} from '../../../modules/local/api/'
include { download_sra }                            from '../../../modules/local/sratools/download_sra'
include { download_tsa }                            from '../../../modules/local/sratools/download_tsa'
include { trinity }                                 from '../../../modules/local/trinity/trinity'
include { gather_files as gather_transcript_files } from '../../../modules/local/gather_files'
include { gather_files as gather_protein_files }    from '../../../modules/local/gather_files'
include { cd_hit_est }                              from '../../../modules/local/cd_hit/run_cd_hit_est'
include { cd_hit as cd_hit1 }                       from '../../../modules/local/cd_hit/run_cd_hit'
include { cd_hit as cd_hit2 }                       from '../../../modules/local/cd_hit/run_cd_hit'
include { download_proteins_in_proteomes }          from '../../../modules/local/uniprot/download_proteins_in_proteomes'
include { download_proteins as download_proteins1 } from '../../../modules/local/uniprot/download_proteins'
include { download_proteins as download_proteins2 } from '../../../modules/local/uniprot/download_proteins'

workflow download_datasets {
    take:
        sra_to_download
        large_protein_set
        close_protein_set
        very_close_protein_set
        transcriptome_set

    main:
        // Transcripts set
        DOWNLOAD_TSA(transcriptome_set)

        DOWNLOAD_SRA(sra_to_download)
        TRINITY(DOWNLOAD_SRA.out)

        DOWNLOAD_TSA.out.join(TRINITY.out)
            .map { id, tsa, trinity -> [id, [tsa[2], trinity[3]], "transcripts_set.fa", 'no']}
            .set { transcriptome_set }
        GATHER_TRANSCRIPT_FILES(transcriptome_set)

        CD_HIT_EST(GATHER_TRANSCRIPT_FILES.out)

        // // Main protein set
        // DOWNLOAD_PROTEINS_IN_PROTEOMES(large_protein_set)
        // DOWNLOAD_PROTEINS1(close_protein_set)

        // DOWNLOAD_PROTEINS_IN_PROTEOMES.out.join(DOWNLOAD_PROTEINS1.out)
        //     .map { id, large, close -> [id, [large[2], close[2]], "main_proteins_set.fa", 'no']}
        //     .set { main_protein_set }
        // GATHER_PROTEIN_FILES(main_protein_set)

        // CD_HIT1(GATHER_PROTEIN_FILES.out.map{id, input -> [id, input, "${id}_main_proteins_set.fa"]})

        // // Parallele protein set
        // DOWNLOAD_PROTEINS2(very_close_protein_set)
        // CD_HIT2(DOWNLOAD_PROTEINS2.out.map{id, input -> [id, input, "${id}_parallele_proteins_set.fa"]})

    // emit:
    transcripts = CD_HIT_EST.out
    // main_proteins = CD_HIT1.out
    // parallele_proteins = CD_HIT2.out
}
