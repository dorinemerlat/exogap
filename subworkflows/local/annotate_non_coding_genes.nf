include { INFERNAL                  } from '../../modules/local/infernal'
include { BARRNAP as BARRNAP_NUCL   } from '../../modules/local/barrnap'
include { BARRNAP as BARRNAP_MITO   } from '../../modules/local/barrnap'
include { RNAMMER                   } from '../../modules/local/rnammer'
// include { AGAT_MERGE_ANNOTATIONS    } from '../../modules/local/agat_merge_annotations'
// include { SPLIT_BY_NCRNA_TYPE       } from '../../modules/local/split_by_ncrna_type'
include { DOWNLOAD_RFAM             } from '../../modules/local/download_rfam'

workflow ANNOTATE_NON_CODING_GENES {
    take:
        genomes

    main:
    DOWNLOAD_RFAM()

    INFERNAL(genomes.combine(DOWNLOAD_RFAM.out.models))
    BARRNAP_NUCL(genomes.map { id, meta, genome -> [id, meta, genome, 'nucl']})
    BARRNAP_MITO(genomes.map { id, meta, genome -> [id, meta, genome, 'mito']})
    RNAMMER(genomes)

    emit:
        infernal = INFERNAL.out
        barrnap_nucl = BARRNAP_NUCL.out
        barrnap_mito = BARRNAP_MITO.out
        rnammer = RNAMMER.out
}
