/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GET TAXONOMIC INFORMATIONS AND SELECT PUBLIC DATA TO USE FOR ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include modules
include { DOWNLOAD_LINEAGE }                from '../../../modules/local/download_lineage.nf'
include { DOWNLOAD_NEWICK }                 from '../../../modules/local/download_newick'
include { SEARCH_PROTEINS_IN_PROTEOMES }    from '../../../modules/local/search_proteins_in_proteomes.nf'
include { SEARCH_PROTEINS }                 from '../../../modules/local/search_proteins'
include { SEARCH_TSA }                      from '../../../modules/local/search_tsa'
include { SEARCH_SRA }                      from '../../../modules/local/search_sra'
include { DOWNLOAD_BUSCO_DATASETS }         from '../../../modules/local/download_busco_datasets'

def gather_genomes(genomes) {
    return genomes
        .map{it -> [it[0], it[1], it[2]]} // remove extra fiels
        .map { id, meta, file -> [id, ["${meta.taxid}": meta.name], file] } // keep only the IDs and essential meta data
        .flatten() // flat and list to get a destructed channel
        .toList()
        .map { genome -> genome.withIndex().collect { it, index -> [index % 3, it] } } // group by two: all the IDs, and all the meta data
        .flatMap { it } // have an it for IDs and one for meta data
        .groupTuple() // group to
        .map { index, it -> it.unique().sort() }
        .toList()
}

def index_genomes_by_lineage(genomes) {
    return genomes
        .map { id, meta, fasta -> [meta.lineage.taxid, id, meta, fasta ] }
        .transpose()
}


def set_selection(genomes, datasets, max_sequence_number) {
    genomes_index_by_lineage = index_genomes_by_lineage(genomes)

    dataset_index_by_taxid = datasets.map { taxid, name, id_list, count, other  -> [taxid, name, id_list, count.readLines()[0].split(',')[2].toInteger(), other ] }
        .filter { it[3] < max_sequence_number } // keep only those a minimal number of proteins
        .filter { it[3] != 0 }

    genomes = genomes_index_by_lineage.combine(dataset_index_by_taxid, by : 0)
        .map { clade_taxid, id, meta, fasta, clade_name, clade_ids, clade_count, clade_other -> [ [id, meta, fasta], [clade_taxid, clade_name, clade_ids, clade_count, clade_other] ]}
        .groupTuple()
        .map { genome, dataset -> [genome, dataset.max { a,b -> a[3] <=> b[3] } ]}
        .map { genome, dataset -> [ genome[0], genome[1], genome[2], ['name': dataset[1], 'taxid': dataset[0], 'id_list': dataset[2], 'count': dataset[3], 'other': dataset[4] ]]}

    return genomes
}

def set_extraction(meta) {
    return meta.map { it -> [it.name, it.taxid, it.id_list] }.unique()
}

workflow GET_INFORMATIONS_ABOUT_GENOMES {
    take:
        genomes

    main:
        // Get the lineage taxonomy of each genome (one by genome)
        DOWNLOAD_LINEAGE(genomes)

        DOWNLOAD_LINEAGE.out
            .map { it -> [ it[0], file(it[1]).splitText() ] }
            .map{ id, lineage -> [ id, lineage.collect { it.split(',') } ] }
            .map{ id, lineage -> [ id, ["lineage": lineage.collect { [ "rank": "${it[0]}", "name": "${it[2]}".replace('\n', ''), "taxid": "${it[1]}" ] }]]}
            .join(genomes)
            .map { id, taxonomy, meta, fasta -> [id, meta + taxonomy, fasta ] }
            .set { genomes }

        DOWNLOAD_NEWICK(gather_genomes(genomes))

        // choose which busco dataset to use for each specie
        DOWNLOAD_BUSCO_DATASETS()

        genomes.map { id, meta, fasta -> meta.lineage.collect{[it.name.toLowerCase(), id, meta, fasta]} } // create an it for each parent of each genome
            .flatMap{ it }
            .combine(DOWNLOAD_BUSCO_DATASETS.out.splitCsv(), by: 0) // -> [taxon, id, meta, fasta, busco_dataset]
            .map { taxon, id, meta, fasta, busco_dataset -> [id, meta, fasta, busco_dataset] }
            .groupTuple()
            .map { id, meta, file, busco_datasets -> [id, meta[0] + [ "busco_dataset": busco_datasets ], file[0]] }
            .set { genomes }

        // select set about data to download (proteins and transcripts)
        genomes.map { id, meta, fasta -> meta.lineage.collect{["taxid": it.taxid, "name": it.name, "rank": it.rank]} } // create an it for each parent of each genome
            .map { it.findAll { it.rank != "superkingdom"} .findAll {it.rank != "species"} } // keep only those with rank different of superkingdom or species
            .flatMap { it }
            .map { [it.taxid, it.name] }
            .unique()
            .set { parents }

        genomes_index_by_lineage = index_genomes_by_lineage(genomes)

        // main proteins set
        SEARCH_PROTEINS_IN_PROTEOMES(parents)
        genomes = set_selection(genomes, SEARCH_PROTEINS_IN_PROTEOMES.out, params.max_proteins_from_a_large_set_of_species)
        genomes.map { genome, meta, fasta, set -> [ genome, meta + [ 'main_proteins_set' : [set.taxid, set] ], fasta]}
            .set { genomes }

        SEARCH_PROTEINS(parents)
        genomes = set_selection(genomes, SEARCH_PROTEINS.out, params.max_proteins_from_relatively_close_species)
        genomes.map { genome, meta, fasta, set -> [ genome, meta + ['main_proteins_set': ['id': [meta.main_proteins_set[0], set.taxid].sort(), 'from_proteomes': meta.main_proteins_set[1], 'from_proteins' : set]], fasta ]}
            .set { genomes }

        // training proteins set
        genomes = set_selection(genomes, SEARCH_PROTEINS.out, params.max_proteins_from_close_species)
        genomes.map { genome, meta, fasta, set -> [ genome, meta + [ 'training_proteins_set' : set ], fasta]}
            .set { genomes }

        // transcripts set
        parents.combine( parents.count() ).set { parents_with_count }
        SEARCH_TSA(parents_with_count)
        genomes = set_selection(genomes, SEARCH_TSA.out, params.max_transcriptomes)
        genomes.map { genome, meta, fasta, set -> [ genome, meta + [ 'transcriptomes_set' : set ], fasta]}
            .set { genomes }

        // extract datasets
        large_proteins_set = set_extraction(genomes.map { id, meta, fasta -> meta.main_proteins_set.from_proteomes })
        small_proteins_set = set_extraction(genomes.map { id, meta, fasta -> meta.main_proteins_set.from_proteins })
        training_proteins_set = set_extraction(genomes.map { id, meta, fasta -> meta.training_proteins_set })
        transcriptome_set = set_extraction(genomes.map { id, meta, fasta -> meta.transcriptomes_set })

        // select complementary transcriptomes (SRA)
        genomes.map { id, meta, fasta -> [[meta.transcriptomes_set.taxid, meta.transcriptomes_set.other.readLines()], [id, meta.taxid]] }
            .filter { !(it[0][1].contains(it[1][1])) }
            .map { tsa_set, genome -> [tsa_set[0], genome[0], genome[1]] }
            .set { sra_to_download }

        sra_to_download.combine( sra_to_download.count() ).set { sra_to_download }

        SEARCH_SRA(sra_to_download)

        SEARCH_SRA.out
            // .map{ clade, specie_name, specie_taxid, out -> [clade, specie_name, specie_taxid, out] }
            .filter { it[3].readLines().size() != 0 }
            .map{ clade, specie_name, specie_taxid, out -> [clade, specie_name, specie_taxid, out, out.readLines().size()] }
            .set { sra_to_download }

    emit:
        genomes                 = genomes

        sra_to_download         = sra_to_download
        transcriptome_set       = transcriptome_set

        large_proteins_set      = large_proteins_set
        small_proteins_set      = small_proteins_set

        training_proteins_set   = training_proteins_set
}
