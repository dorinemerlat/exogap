include { GENERATE_RUNLIST } from '../../../modules/local/module_genes/generate_runlist.nf'
include { VARUS } from '../../../modules/local/module_genes/varus.nf'
include { BRAKER } from '../../../modules/local/module_genes/braker.nf'
include { GALBA } from '../../../modules/local/module_genes/galba.nf'

workflow GENES_ANNOTATION {
    take:
        masked_genomes
    main:

        // filter masked_genomes to keep only those where meta.SRA is not null
        masked_genomes
            .map { id, meta, fasta -> [id, meta, meta.lineage.species.name] }
            .set { genomes_with_sra }

        GENERATE_RUNLIST(genomes_with_sra)

        unmasked_genomes
            .combine(GENERATE_RUNLIST.out, by: 0)
            .map { id, meta1, fasta, meta2, runlist -> [id, meta1, fasta, runlist, meta1.lineage.species.name] }
            .set { varus_input }

        // VARUS(varus_input)

        masked_genomes
            .join(VARUS.out.bam, by: 0)
            .map { id, meta1, fasta, meta2, bam -> [id, meta1, fasta, bam, file(params.protein_set), meta1.lineage.species.name] }
            .branch { id, meta, fasta, bam, proteins, specie_name ->
                with_bam: bam.size() > 0
                without_bam: bam.size() == 0
            }
            .set { annotation_input }

        annotation_input.with_bam
            .set { braker_input }

        annotation_input.without_bam
            .map { id, meta, fasta, bam, proteins, specie_name -> [id, meta, fasta, proteins, specie_name] }
            .set { galba_input }

        // annotation with BRAKER if RNA-seq data is available, otherwise with GALBA
        BRAKER(annotation_input.with_bam)
        // GALBA(galba_input)
    // emit:
}

