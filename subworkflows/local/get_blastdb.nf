/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GET TAXONOMIC INFORMATIONS AND SELECT PUBLIC DATA TO USE FOR ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include modules
include { UPDATE_LOCALDB }                    from '../../modules/local/update_localdb'

workflow GET_BLASTDB {
    take:
        blastdb

    main:
        UPDATE_LOCALDB(blastdb)

        UPDATE_LOCALDB.map { it -> [it.find { it =~ ".phr" }.toString().replaceFirst(/.phr/, ""), it] }
            .set { blastdb }

    emit:
        blastdb = blastdb

}
