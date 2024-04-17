include { AGAT_MERGE_ANNOTATIONS    } from '../../modules/local/agat_merge_annotations'

workflow POSTPROCESS {
    take:
        genomes
        protein_coding_genes_gff // id, meta, gff
        infernal_gff // id, meta, gff
        barrnap_nucl_gff
        barrnap_mito_gff
        rnammer_gff
        repeats

    main:
    genomes
        .join(protein_coding_genes_gff, by: [0,1])
        .join(infernal_gff, by: [0,1])
        .join(barrnap_nucl_gff, by: [0,1])
        .join(barrnap_mito_gff, by: [0,1])
        .join(rnammer_gff, by: [0,1])
        .join(repeats, by: [0,1])
        .map { id, meta, genomes, proteins, infernal, barrnap_nucl, barrnap_mito, rnammer, repeats -> [id, meta, genomes, [proteins, infernal, barrnap_nucl, barrnap_mito, rnammer, repeats], 'postprocess']}
        .set { genomes_for_merging }

    AGAT_MERGE_ANNOTATIONS(genomes_for_merging)

    emit:
        gff = AGAT_MERGE_ANNOTATIONS.out
}
