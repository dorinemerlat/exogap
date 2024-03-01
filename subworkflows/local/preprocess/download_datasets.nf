/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DOWNLOAD SELECTED PUBLIC DATA TO USE FOR ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include modules

// include {} from '../../../modules/local/api/'

workflow DOWNLOAD_DATASETS {
    take:
        genomes
        sra_to_download
        large_protein_set
        close_protein_set
        very_close_protein_set
        transcriptome_set

    main:
        genomes.view()
        // download the selected proteins (large set, close set, very close)
        // DOWNLOAD_PROTEINS()
        // download the selected transcriptomes

        // download the selected SRA data

        // add SRA to transcriptomes

        // eliminate redundancy in transcripts sets

        // eliminate redundancy in main proteins sets
        // 1. Concatenate large protein set + close protein set
        // 2. Eliminate redundancy in concatenated set

        // eliminate redundancy in closed proteins set


    // emit:

}
