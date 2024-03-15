/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GET TAXONOMIC INFORMATIONS AND SELECT PUBLIC DATA TO USE FOR ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include modules
include { DOWNLOAD_LINEAGE }                                from '../../../modules/local/download_lineage.nf'
include { DOWNLOAD_NEWICK }                                 from '../../../modules/local/download_newick'
include { SEARCH_PROTEINS_IN_PROTEOMES }                    from '../../../modules/local/search_proteins_in_proteomes.nf'
include { SEARCH_PROTEINS }                                 from '../../../modules/local/search_proteins'
include { SEARCH_TSA }                                      from '../../../modules/local/search_tsa'
include { SEARCH_SRA }                                      from '../../../modules/local/search_sra'
include { DOWNLOAD_BUSCO_DATASETS }                         from '../../../modules/local/download_busco_datasets'
include { DOWNLOAD_SRA }                                    from '../../../modules/local/download_sra'
include { DOWNLOAD_TSA }                                    from '../../../modules/local/download_tsa'
include { TRINITY }                                         from '../../../modules/local/trinity'
include { GATHER_FILES as GATHER_TRANSCRIPT_FILES }         from '../../../modules/local/gather_files'
include { GATHER_FILES as GATHER_MAIN_PROTEIN_FILES }       from '../../../modules/local/gather_files'
include { GATHER_FILES as GATHER_TRAINING_PROTEIN_FILES }   from '../../../modules/local/gather_files'
include { CD_HIT_EST }                                      from '../../../modules/local/cd_hit_est'
include { CD_HIT as CD_HIT_MAIN_PROTEINS }                  from '../../../modules/local/cd_hit'
include { CD_HIT as CD_HIT_TRAINING_PROTEINS }              from '../../../modules/local/cd_hit'
include { DOWNLOAD_PROTEINS_IN_PROTEOMES }                  from '../../../modules/local/download_proteins_in_proteomes'
include { DOWNLOAD_PROTEINS as DOWNLOAD_PROTEINS1 }         from '../../../modules/local/download_proteins'
include { DOWNLOAD_PROTEINS as DOWNLOAD_PROTEINS2 }         from '../../../modules/local/download_proteins'

def setSelection(genomes, datasets, max_sequence_number) {
    genomes_index_by_lineage = Utils.indexGenomesByLineage(genomes)

    dataset_index_by_taxid = datasets.map { taxid, name, id_list, count, other  -> [taxid, name, id_list, count.readLines()[0].split(',')[2].toInteger(), other ] }
        .filter { it[3] < max_sequence_number } // keep only those a minimal number of proteins
        .filter { it[3] != 0 }


    genomes = genomes_index_by_lineage.combine(dataset_index_by_taxid, by : 0)
        .map { clade_taxid, id, meta, fasta, clade_name, clade_ids, clade_count, clade_other -> [ [id, meta, fasta], [clade_taxid, clade_name, clade_ids, clade_count, clade_other] ]}
        .groupTuple()
        .map { genome, dataset -> [genome, dataset.max { a,b -> a[3] <=> b[3] } ]}
        .map { genome, dataset -> [ dataset[0], ['taxid': dataset[0], 'name': dataset[1], 'count': dataset[3], 'other': dataset[4]], dataset[2]]}
        .unique()

    return genomes
}

def formatName(name) {
    // remove '[', ']', ',')'
    return name.join('_').replaceAll('[\\[\\],]', '').replaceAll(' ', '_')
}
def getSetFiles(genomes, output_name, header) {
    return genomes.groupTuple()
        .map { set, genome_id -> [formatName(set[1]),  ['name': formatName(set[1]), 'taxid': formatName(set[0]), 'genomes': genome_id], set[2].flatten(), output_name, header] }
}

def concatenateGenomeAndSet(genomes, dataset){
    return dataset.map { name, meta, file -> [meta.genomes, [name, meta.taxid, file], meta.genomes.size()]}
                .map { genome, set, size -> [genome, Collections.nCopies(size, set)] }
                .transpose()
                .join( genomes )
                .map { id, set, genome, fasta -> [id, genome, fasta, ['name':set[0], 'taxid': set[1], 'file': set[2]]] }
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

        DOWNLOAD_NEWICK(Utils.gatherGenomes(genomes))

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

        // select main proteins set
        if (params.max_proteins_for_large_main != null && params.max_proteins_for_large_main != false) {
            SEARCH_PROTEINS_IN_PROTEOMES(parents)
            DOWNLOAD_PROTEINS_IN_PROTEOMES(setSelection(genomes, SEARCH_PROTEINS_IN_PROTEOMES.out, params.max_proteins_for_large_main))

            Utils.indexGenomesByLineage(genomes).combine(DOWNLOAD_PROTEINS_IN_PROTEOMES.out, by : 0)
                .map { set_taxid, genome_id, genome_meta, genome_fasta, set_meta, set_file -> [genome_id, Utils.updateLinkedHashMap(
                                                                                                    genome_meta, 'main_protein_set', Utils.updateLinkedHashMap(
                                                                                                    genome_meta.main_protein_set, 'from_uniprot_proteomes',[set_taxid, set_meta, set_file] )), genome_fasta]}
                .set { genomes }


        } else {
            genomes.map { genome, meta, fasta, set -> [ genome, Utils.updateLinkedHashMap(meta, 'main_protein_set', Utils.updateLinkedHashMap(meta.main_protein_set, 'from_uniprot_proteomes', Utils.createEmptySet())), fasta]}
                .set { genomes }
        }

        if (params.max_proteins_for_close_main != null && params.max_proteins_for_close_main != false) {
            SEARCH_PROTEINS(parents)
            DOWNLOAD_PROTEINS1(setSelection(genomes, SEARCH_PROTEINS.out, params.max_proteins_for_close_main))

            Utils.indexGenomesByLineage(genomes).combine(DOWNLOAD_PROTEINS1.out, by : 0)
                .map { set_taxid, genome_id, genome_meta, genome_fasta, set_meta, set_file -> [genome_id, Utils.updateLinkedHashMap(
                                                                                                    genome_meta, 'main_protein_set', Utils.updateLinkedHashMap(
                                                                                                    genome_meta.main_protein_set, 'from_uniprot_proteins',[set_taxid, set_meta, set_file] )), genome_fasta]}
                .set { genomes }

        } else {
            genomes.map { genome, meta, fasta, set -> [ genome, Utils.updateLinkedHashMap(meta, 'main_protein_set', Utils.updateLinkedHashMap(meta.main_protein_set, 'from_uniprot_proteins', Utils.createEmptySet())), fasta]}
                .set { genomes }
        }

        // concatenate downloaded proteins from proteomes, download proteins and personal proteins
        main_proteins = getSetFiles(genomes.map { id, meta, fasta -> [id, meta, meta.main_protein_set.from_uniprot_proteomes, meta.main_protein_set.from_uniprot_proteins, meta.main_protein_set.personal_set] }
            .map { id, meta, from_uniprot_proteomes, from_uniprot_proteins, personal_set -> [[
                                                                [from_uniprot_proteomes[0], from_uniprot_proteins[0], Utils.getBaseNameOrNull(personal_set)].findAll { it != null }, // id
                                                                [from_uniprot_proteomes[1].name, from_uniprot_proteins[1].name, Utils.getBaseNameOrNull(personal_set)].findAll { it != null }, // name
                                                                [from_uniprot_proteomes[2], from_uniprot_proteins[2], personal_set].findAll { it != null }], // file
                                                                id] }, "main_proteins_set.fa", 'no')

    GATHER_MAIN_PROTEIN_FILES(main_proteins)
        CD_HIT_MAIN_PROTEINS(GATHER_MAIN_PROTEIN_FILES.out)
        genomes = concatenateGenomeAndSet(genomes, CD_HIT_MAIN_PROTEINS.out)
        genomes.map { genome, meta, fasta, set -> [ genome, Utils.updateLinkedHashMap(meta, 'main_protein_set', Utils.updateLinkedHashMap(meta.main_protein_set, 'dataset', set)), fasta]}
                .set { genomes }

        // training proteins set
        if (params.max_proteins_for_training != null && params.max_proteins_for_training != false) {
            DOWNLOAD_PROTEINS2(setSelection(genomes, SEARCH_PROTEINS.out, params.max_proteins_for_training))

            Utils.indexGenomesByLineage(genomes).combine(DOWNLOAD_PROTEINS2.out, by : 0)
                .map { set_taxid, genome_id, genome_meta, genome_fasta, set_meta, set_file -> [genome_id, Utils.updateLinkedHashMap(
                                                                                                    genome_meta, 'training_protein_set', Utils.updateLinkedHashMap(
                                                                                                    genome_meta.training_protein_set, 'from_uniprot',[set_taxid, set_meta, set_file] )), genome_fasta]}
                .set { genomes }

        } else {
            genomes.map { genome, meta, fasta, set -> [ genome, Utils.updateLinkedHashMap(meta, 'training_protein_set', Utils.updateLinkedHashMap(meta.training_protein_set, 'from_uniprot', Utils.createEmptySet())), fasta]}
                .set { genomes }
        }

        // concatenate download proteins and personal proteins
        training_proteins = getSetFiles(genomes.map { id, meta, fasta -> [id, meta, meta.training_protein_set.from_uniprot, meta.training_protein_set.personal_set] }
            .map { id, meta, from_uniprot, personal_set -> [[
                                                                [from_uniprot[0], Utils.getBaseNameOrNull(personal_set)].findAll { it != null }, // id
                                                                [from_uniprot[1].name, Utils.getBaseNameOrNull(personal_set)].findAll { it != null }, // name
                                                                [from_uniprot[2], personal_set].findAll { it != null }], // file
                                                                id] }, "training_protein_set.fa", 'no')

        GATHER_TRAINING_PROTEIN_FILES(training_proteins)
        CD_HIT_TRAINING_PROTEINS(GATHER_TRAINING_PROTEIN_FILES.out)
        genomes = concatenateGenomeAndSet(genomes, CD_HIT_TRAINING_PROTEINS.out)
        genomes.map { genome, meta, fasta, set -> [ genome, Utils.updateLinkedHashMap(meta, 'training_protein_set', Utils.updateLinkedHashMap(meta.training_protein_set, 'dataset', set)), fasta]}
                .set { genomes }

        // transcripts set
        if (params.max_transcriptomes != null && params.max_transcriptomes != false) {
            parents.combine( parents.count() ).set { parents_with_count }
            SEARCH_TSA(parents_with_count)
            DOWNLOAD_TSA(setSelection(genomes, SEARCH_TSA.out, params.max_transcriptomes))

            Utils.indexGenomesByLineage(genomes).combine(DOWNLOAD_TSA.out, by : 0)
                .map { set_taxid, genome_id, genome_meta, genome_fasta, set_meta, set_file -> [genome_id, Utils.updateLinkedHashMap(
                                                                                                    genome_meta, 'transcript_set', Utils.updateLinkedHashMap(
                                                                                                    genome_meta.transcript_set, 'from_tsa',[set_taxid, set_meta, set_file] )), genome_fasta]}
                .set { genomes }

        } else {
            genomes.map { genome, meta, fasta, set -> [ genome, Utils.updateLinkedHashMap(meta, 'transcript_set', Utils.updateLinkedHashMap(meta.transcript_set, 'from_tsa', Utils.createEmptySet())), fasta]}
                .set { genomes }
        }


        // select complementary transcriptomes (SRA)
        if (params.use_sra == true) {
            genomes.map { id, meta, fasta -> [[meta.transcript_set.from_tsa[1].taxid, meta.transcript_set.from_tsa[1].other.readLines()], [id, meta.taxid]] }
                .filter { !(it[0][1].contains(it[1][1])) }
                .map { tsa_set, genome -> [tsa_set[0], genome[0], genome[1]] }
                .set { sra_to_download }

            sra_to_download.combine( sra_to_download.count() ).set { sra_to_download }
            SEARCH_SRA(sra_to_download)

            SEARCH_SRA.out
                .filter { it[3].readLines().size() != 0 }
                .map{ clade, specie_name, specie_taxid, out -> [clade, specie_name, specie_taxid, out, out.readLines().size()] }
                .set { sra_to_download }

                DOWNLOAD_SRA(sra_to_download)
                TRINITY(DOWNLOAD_SRA.out)

                Utils.indexGenomesByLineage(genomes).combine(TRINITY.out.groupTuple().map {id, name, taxid, file -> [id, ['taxid': taxid, 'name': name, 'count': taxid.size()], file] }, by : 0)
                .map { set_taxid, genome_id, genome_meta, genome_fasta, set_meta, set_file -> [genome_id, Utils.updateLinkedHashMap(
                                                                                                    genome_meta, 'transcript_set', Utils.updateLinkedHashMap(
                                                                                                    genome_meta.transcript_set, 'from_sra',[set_taxid, set_meta, set_file] )), genome_fasta]}
                .set { genomes }

        } else {
            genomes.map { genome, meta, fasta, set -> [ genome, Utils.updateLinkedHashMap(meta, 'transcript_set', Utils.updateLinkedHashMap(meta.transcript_set, 'from_sra', Utils.createEmptySet())), fasta]}
                .set { genomes }
        }

        // concatenate download tsa, download_sra and personal transcripts
        transcripts = getSetFiles(genomes.map { id, meta, fasta -> [id, meta, meta.transcript_set.from_tsa, meta.transcript_set.from_sra, meta.transcript_set.personal_set] }
            .map { id, meta, from_tsa, from_sra, personal_set -> [[
                                                                [from_tsa[0], from_sra[0], Utils.getBaseNameOrNull(personal_set)].findAll { it != null }, // id
                                                                [from_tsa[1].name, from_sra[1].name, Utils.getBaseNameOrNull(personal_set)].findAll { it != null }, // name
                                                                [from_tsa[2], from_sra[2], personal_set].findAll { it != null }], // file
                                                                id] }, "transcript_set.fa", 'no')
        GATHER_TRANSCRIPT_FILES(transcripts)
        CD_HIT_EST(GATHER_TRANSCRIPT_FILES.out.map { id, meta, file -> [id, meta, file, '5', '10', '0'] })
        genomes = concatenateGenomeAndSet(genomes, CD_HIT_EST.out)
        genomes.map { genome, meta, fasta, set -> [ genome, Utils.updateLinkedHashMap(meta, 'transcript_set', Utils.updateLinkedHashMap(meta.transcript_set, 'dataset', set)), fasta]}
                .set { genomes }

    emit:
        genomes = genomes
}
