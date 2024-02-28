/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GET TAXONOMIC INFORMATIONS AND SELECT PUBLIC DATA TO USE FOR ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include modules
include { GET_TAXONOMIC_LINEAGE }       from '../../../modules/local/api/get_taxonomic_lineage'
include { GET_NEWICK }                  from '../../../modules/local/api/get_newick'
include { COUNT_PROTEINS_IN_PROTEOMES } from '../../../modules/local/api/count_proteins_in_proteomes'
include { COUNT_PROTEINS }              from '../../../modules/local/api/count_proteins'
include { COUNT_TRANSCRIPTOMES }        from '../../../modules/local/api/count_transcriptomes'
include { GATHER_PUBLIC_DATA }          from '../../../modules/local/gather_public_data'
include { DOWNLOAD_BUSCO_DATASETS }     from '../../../modules/local/busco/download_busco_datasets'
include { SELECT_PUBLIC_DATA }          from '../../../modules/local/select_public_data'

def gather_genomes(genomes) {
    return genomes
        .map{it -> [it[0], it[1], it[2]]} // remove extra fiels
        .map { id, meta, file -> [id, ["${meta.taxid}": meta.name], file] } // keep only the IDs and essential meta data
        .flatten() // flat and list to get a destructed channel
        .toList()
        .map { genome -> genome.withIndex().collect { it, index -> [index % 3, it] } } // group by two: all the IDs, and all the meta data
        .flatMap { it } // have an item for IDs and one for meta data
        .groupTuple() // group to
        .map { index, it -> it.unique().sort() }
        .toList()
}

workflow GET_INFORMATIONS_ABOUT_GENOMES {
    take:
        genomes

    main:
        // Get the lineage taxonomy of each genome (one by genome)
        GET_TAXONOMIC_LINEAGE(genomes)

        GET_TAXONOMIC_LINEAGE.out
            .map { it -> [ it[0], file(it[1]).splitText() ] }
            .map{ id, lineage -> [ id, lineage.collect { it.split(',') } ] }
            .map{ id, lineage -> [ id, ["lineage": lineage.collect { [ "rank": "${it[0]}", "name": "${it[2]}".replace('\n', ''), "taxid": "${it[1]}" ] }]]}
            .join(genomes)
            .map { id, taxonomy, meta, fasta -> [id, meta + taxonomy, fasta ] }
            .set { genomes }

        GET_NEWICK(gather_genomes(genomes))

        // choose which busco dataset to use for each specie
        DOWNLOAD_BUSCO_DATASETS()

        genomes.map { id, meta, fasta -> meta.lineage.collect{[it.name.toLowerCase(), id, meta, fasta]} } // create an item for each parent of each genome
            .flatMap{ it }
            .combine(DOWNLOAD_BUSCO_DATASETS.out.splitCsv(), by: 0) // -> [taxon, id, meta, fasta, busco_dataset]
            .map { taxon, id, meta, fasta, busco_dataset -> [id, meta, fasta, busco_dataset] }
            .groupTuple()
            .map { id, meta, file, busco_datasets -> [id, meta[0] + [ "busco_dataset": busco_datasets ], file[0]] }
            .set { genomes }

        // select set about data to download
        genomes.map { id, meta, fasta -> meta.lineage.collect{["taxid": it.taxid, "name": it.name, "rank": it.rank]} } // create an item for each parent of each genome
            .map { it.findAll { it.rank != "superkingdom"} .findAll {it.rank != "species"} } // keep only those with rank different of superkingdom or species
            .flatMap { it }
            .map { [it.taxid, it.name] }
            .unique()
            .set { parents }

        COUNT_PROTEINS_IN_PROTEOMES(parents)

        COUNT_PROTEINS_IN_PROTEOMES.out.proteins_count
            .map { taxid, name, count -> [taxid, name, count.readLines()[0].split(',')[2].toInteger() ] }
            .filter { it[2] < params.max_proteins_from_a_large_set_of_species } // keep only those a minimal number of proteins
            .max { it[2] }
            .view()

        COUNT_PROTEINS(parents)

        COUNT_PROTEINS.out.proteins_count
            .map { taxid, name, count -> [taxid, name, count.readLines()[0].split(',')[2].toInteger()] }
            .filter { it[2] < params.max_proteins_from_relatively_close_species } // keep only those a minimal number of proteins
            .max { it[2] }
            .view()

        parents
            .combine( parents.count() )
            .set { parents_with_count }

        COUNT_TRANSCRIPTOMES(parents_with_count)

        COUNT_TRANSCRIPTOMES.out.transcriptomes_count
            .map { taxid, name, count -> [taxid, name, count.readLines()[0].split(',')[2].toInteger() ] }
            .view()

        // COUNT_PROTEINS_IN_PROTEOMES.out.view()
        // GATHER_PUBLIC_DATA(
        //     COUNT_PROTEINS_PROTEOMES.out.proteins_count.map {taxid, name, csv -> csv }.collect(),
        //     COUNT_PROTEINS.out.proteins_count.map {taxid, name, csv -> csv }.collect(),
        //     COUNT_TRANSCRIPTOMES.out.transcriptomes_count.map {taxid, name, csv -> csv }.collect())

        // genomes
        //     .map {id, meta, fasta -> [ id, meta, fasta, meta.lineage.collect{it. taxid}.join(" ") ] }
        //     .combine (GATHER_PUBLIC_DATA.out)
        //     .set { genomes_with_taxonomic_informations }

        // SELECT_PUBLIC_DATA(genomes_with_taxonomic_informations)

    emit:
        genomes                   = genomes
    //     proteins_from_proteomes   = SELECT_PROTEINS_FROM_PROTEOMES.out
    //     proteins_from_closest     = SELECT_PROTEINS_FROM_CLOSEST_SPECIES.out
    //     transcripts_from_sra      = SELECT_TRANSCRIPTS_FROM_SRA.out
    //     transcripts_from_tsa      = SELECT_TRANSCRIPTS_FROM_TSA.out

}
