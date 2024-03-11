/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DOWNLOAD SELECTED PUBLIC DATA TO USE FOR ANNOTATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include modules

// include {} from '../../../modules/local/api/'
include { DOWNLOAD_SRA }                            from '../../../modules/local/download_sra'
include { DOWNLOAD_TSA }                            from '../../../modules/local/download_tsa'
include { TRINITY }                                 from '../../../modules/local/trinity'
include { GATHER_FILES as GATHER_TRANSCRIPT_FILES } from '../../../modules/local/gather_files'
include { GATHER_FILES as GATHER_PROTEIN_FILES }    from '../../../modules/local/gather_files'
include { CD_HIT_EST }                              from '../../../modules/local/cd_hit_est'
include { CD_HIT as CD_HIT1 }                       from '../../../modules/local/cd_hit'
include { CD_HIT as CD_HIT2 }                       from '../../../modules/local/cd_hit'
include { DOWNLOAD_PROTEINS_IN_PROTEOMES }          from '../../../modules/local/download_proteins_in_proteomes'
include { DOWNLOAD_PROTEINS as DOWNLOAD_PROTEINS1 } from '../../../modules/local/download_proteins'
include { DOWNLOAD_PROTEINS as DOWNLOAD_PROTEINS2 } from '../../../modules/local/download_proteins'


//  def set_selection(genomes, datasets, dataset_name, max_sequence_number) {
//     genomes_index_by_dataset = genomes.map { id, meta, fasta  -> [ meta.proteins_in_proteomes_set.taxid, id meta, fasta ] }

//     genomes_index_by_dataset.view()

//     // genomes = genomes_index_by_lineage.combine(dataset_index_by_taxid, by : 0)
//     //     .map { clade_taxid, id, meta, fasta, clade_name, clade_ids, clade_count, clade_other -> [ [id, meta, fasta], [clade_taxid, clade_name, clade_ids, clade_count, clade_other] ]}
//     //     .groupTuple()
//     //     .map { genome, dataset -> [genome, dataset.max { a,b -> a[3] <=> b[3] } ]}
//     //     .map { genome, dataset -> [ genome[0], genome[1], genome[2], ['name': dataset[1], 'taxid': dataset[0], 'id_list': dataset[2], 'count': dataset[3], 'other': dataset[4] ]]}

//     // if (dataset_name == "proteins_in_proteomes_set") {
//     //     return genomes.map { genome, meta, fasta, set -> [ genome, meta + [ proteins_in_proteomes_set : set ], fasta]}
//     // } else if (dataset_name == "proteins_set") {
//     //     return genomes.map { genome, meta, fasta, set -> [ genome, meta + [ proteins_set : set ], fasta]}
//     // } else if (dataset_name == "closed_proteins_set") {
//     //     return genomes.map { genome, meta, fasta, set -> [ genome, meta + [ closed_proteins_set : set ], fasta]}
//     // } else if (dataset_name == "transcriptomes_set") {
//     //     return genomes.map { genome, meta, fasta, set -> [ genome, meta + [ transcriptomes_set : set ], fasta]}
//     // }
// }


workflow DOWNLOAD_DATASETS {
    take:
        genomes
        sra_to_download
        transcriptome_set
        large_protein_set
        small_proteins_set
        training_proteins_set

    main:
        // Transcripts set
        DOWNLOAD_TSA(transcriptome_set)

        DOWNLOAD_SRA(sra_to_download)
        TRINITY(DOWNLOAD_SRA.out)

        DOWNLOAD_TSA.out.map { name, id, file -> [id, name, file ]}
            .combine (TRINITY.out, by: 0)
            .groupTuple()
            .map { id, clade, tsa, specie_names, specie_taxids, sras -> [id, tsa.unique() + sras]}
            .map { id, tsa -> [id, tsa, "transcripts_set.fa", 'no']}
            .set { transcriptome_set }

        GATHER_TRANSCRIPT_FILES(transcriptome_set)
        CD_HIT_EST(GATHER_TRANSCRIPT_FILES.out)

        // set path to transcripts set in genomes's meta
        genomes.map { id, meta, fasta  -> [ meta.transcriptomes_set.taxid, [id, meta, fasta] ] }
            .combine ( CD_HIT_EST.out, by: 0 )
            .map { id, genome, dataset_file -> [ genome[0], genome[1] + ['transcriptomes_set': genome[1].transcriptomes_set + [ 'dataset' : dataset_file]], genome[2]] }
            .set { genomes }

        // Main protein set
        genomes.map { id, meta, fasta -> [meta.main_proteins_set.id, [id, meta, fasta]] }
            .set { genomes_index_by_main_proteins_set }

        DOWNLOAD_PROTEINS_IN_PROTEOMES(large_protein_set)
        DOWNLOAD_PROTEINS1(small_proteins_set)

        DOWNLOAD_PROTEINS_IN_PROTEOMES.out.map { id, name, fasta -> [[id, fasta]] }
            .combine ( DOWNLOAD_PROTEINS1.out.map { id, name, fasta -> [[id, fasta]] } )
            .map { large, small -> [[large[0], small[0]].sort(), [large[1], small[1]]] }
            .combine ( genomes_index_by_main_proteins_set, by: 0 )
            .map { id, files, genomes -> [ id, files, "main_proteins_set.fa", 'no' ]}
            .unique()
            .set { main_protein_set_proteins }

        GATHER_PROTEIN_FILES(main_protein_set_proteins)
        CD_HIT1(GATHER_PROTEIN_FILES.out)

        // set path to main proteins set  in genomes's meta
        genomes_index_by_main_proteins_set
            .combine ( CD_HIT1.out, by: 0 )
            .map { id, genome, dataset_file -> [ genome[0], genome[1] + ['main_proteins_set': genome[1].main_proteins_set + [ 'dataset' : dataset_file]], genome[2]] }
            .set { genomes }

        // training protein set
        DOWNLOAD_PROTEINS2(training_proteins_set)
        CD_HIT2(DOWNLOAD_PROTEINS2.out.map { id, name, file -> [id, file] })

        // set path to training proteins set in genomes's meta
        genomes.map { id, meta, fasta  -> [ meta.training_proteins_set.taxid, [id, meta, fasta] ] }
            .combine ( CD_HIT2.out, by: 0 )
            .map { id, genome, dataset_file -> [ genome[0], genome[1] + ['training_proteins_set': genome[1].training_proteins_set + [ 'dataset' : dataset_file]], genome[2]] }
            .set { genomes }


    emit:
        genomes = genomes
}
