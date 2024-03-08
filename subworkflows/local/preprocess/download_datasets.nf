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
        download_tsa(transcriptome_set)

        download_sra(sra_to_download)
        trinity(download_sra.out)

        download_tsa.out.join(trinity.out)
            .map { id, tsa, trinity -> [id, [tsa[2], trinity[3]], "transcripts_set.fa", 'no']}
            .set { transcriptome_set }
        gather_transcript_files(transcriptome_set)

        cd_hit_est(gather_transcript_files.out)

        // // Main protein set
        // download_proteins_in_proteomes(large_protein_set)
        // download_proteins1(close_protein_set)

        // download_proteins_in_proteomes.out.join(download_proteins1.out)
        //     .map { id, large, close -> [id, [large[2], close[2]], "main_proteins_set.fa", 'no']}
        //     .set { main_protein_set }
        // gather_protein_files(main_protein_set)

        // cd_hit1(gather_protein_files.out.map{id, input -> [id, input, "${id}_main_proteins_set.fa"]})

        // // Parallele protein set
        // download_proteins2(very_close_protein_set)
        // cd_hit2(download_proteins2.out.map{id, input -> [id, input, "${id}_parallele_proteins_set.fa"]})

    // emit:
    transcripts = cd_hit_est.out
    // main_proteins = cd_hit1.out
    // parallele_proteins = cd_hit2.out
}
