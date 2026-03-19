include { DOWNLOAD_SRA } from "$projectDir/modules/local/module_genes/download_sra.nf"
include { STAR_SINGLE } from "$projectDir/modules/local/module_genes/star_single.nf"
include { STAR_PAIRED } from "$projectDir/modules/local/module_genes/star_paired.nf"
include { MERGE_BAM } from "$projectDir/modules/local/module_genes/merge_bam.nf"
include { BRAKER2 } from "$projectDir/modules/local/module_genes/braker2.nf"
include { BRAKER3 } from "$projectDir/modules/local/module_genes/braker3.nf"
include { GALBA } from "$projectDir/modules/local/module_genes/galba.nf"
include { AGAT_STATISTICS } from "$projectDir/modules/local/module_genes/agat_statistics.nf"
include { AGAT_FILTER_GENE_BY_LENGTH } from "$projectDir/modules/local/module_genes/agat_filter_gene_by_length.nf"
include { AGAT_FILTER_INCOMPLETE_GENES } from "$projectDir/modules/local/module_genes/agat_filter_incomplete_genes.nf"
include { EXTRACT_CANONICAL_PROTEINS } from "$projectDir/modules/local/module_genes/extract_canonical_proteins.nf"
include { OMARK } from "$projectDir/modules/local/module_genes/omark.nf"
include { BUSCO } from "$projectDir/modules/local/module_genes/busco.nf"
include { FILTER_OMARK_CONTAMINATIONS } from "$projectDir/modules/local/module_genes/filter_omark_contaminations.nf"
include { OMARK as OMARK_ON_CONTAMINATION_FILTERED } from "$projectDir/modules/local/module_genes/omark.nf"
include { BUSCO as BUSCO_ON_CONTAMINATION_FILTERED } from "$projectDir/modules/local/module_genes/busco.nf"
include { CLASSIFY_GENES_VS_REPEATS } from "$projectDir/modules/local/module_genes/classify_genes_vs_repeats.nf"
// include { INTERPROSCAN } from "$projectDir/modules/local/module_genes/interproscan.nf"
// include { BLASTP } from "$projectDir/modules/local/module_genes/blastp.nf"

workflow GENES_ANNOTATION {
    take:
        masked_genomes
        repeats_gff

    main:
        // check if RNASeq-file is provided in the metadata, if not set rnaseq_runs to null
        masked_genomes
            .filter { id, meta, fasta -> meta["rnaseq_runs"] != null }
            .flatMap { id, meta, fasta ->
                meta["rnaseq_runs"].collect { run_tuple ->
                    def run_files = run_tuple.drop(1).findAll { it != null }
                    [id, run_files]
                }
            }
            .set { rnaseq_runs_files }

        // filter masked_genomes to keep only those where meta.SRA is not null
        masked_genomes
            .map { id, meta, fasta -> [id, meta, meta.lineage.species.name] }
            .set { genomes_with_sra }

        GENERATE_RUNLIST(genomes_with_sra)

        // filter out empty runlists
        GENERATE_RUNLIST.out
            .filter { id, meta, file ->
                file.readLines().size() > 1
            }
            .set { runlists_non_empty }

        runlists_non_empty
            .flatMap { id, meta, runlist ->
                runlist.readLines()
                    .drop(1)
                    .findAll { it?.trim() }
                    .collect { line ->
                        def cols = line.trim().split(/\s+/)
                        [id, cols[0]]
                    }
            }
            .set { runlist_sra_ids }

        DOWNLOAD_SRA(runlist_sra_ids)

        DOWNLOAD_SRA.out
            .concat( rnaseq_runs_files )
            .branch { id, sra_files ->
                    def files = (sra_files instanceof List) ? sra_files : [sra_files]

                    single: files.size() == 1
                    paired: files.size() == 2
                }
            .set { star_input }

        star_input.single
            .groupTuple(by: 0)
            .join(masked_genomes, by: 0)
            .map { id, fastq, meta, genome -> [id, meta, genome, fastq.flatten()] }
            .set { star_single_input }

        star_input.paired
            .groupTuple(by: 0)
            .join(masked_genomes, by: 0)
            .map { id, fastq, meta, genome -> [id, meta, genome, fastq.flatten()] }
            .set { star_paired_input }

        STAR_SINGLE(star_single_input)
        STAR_PAIRED(star_paired_input)

        STAR_SINGLE.out.bam
            .concat(STAR_PAIRED.out.bam)
            .groupTuple(by: [0, 1])
            .map { id, meta, bams -> [id, meta, bams ] }
            .branch { id, meta, bams ->
                has_multiple_bams: bams.size() > 1
                single_bam: bams.size() == 1
            }
            .set { merged_star_input }

        MERGE_BAM(merged_star_input.has_multiple_bams)

        MERGE_BAM.out
            .concat( merged_star_input.single_bam )
            .set { all_bams }

        masked_genomes
            .join(all_bams, by: [0, 1])
            .map { id, meta, genome, bam -> [id, meta, "star", genome, bam[0], file(params.protein_set), meta.lineage.species.name] }
            .set { braker_with_star_input }

        masked_genomes
            .join(all_bams, by:[0, 1])
            .map { id, meta, genome, bam -> [id, meta, "varus", genome, bam[0], file(params.protein_set), meta.lineage.species.name] }
            .set { braker_with_varus_input }

        braker_with_star_input
            .concat(braker_with_varus_input)
            .set { braker3_input }

        BRAKER3(braker3_input)

        masked_genomes
            .map { id, meta, genome -> [id, meta, genome, file(params.protein_set), meta.lineage.species.name] }
            .set { braker2_and_galba_input }

        BRAKER2(braker2_and_galba_input)
        GALBA(braker2_and_galba_input)

        BRAKER3.out.gff.map { id, meta, rnaseq_aligner, gff -> [id, meta, "braker_with_" + rnaseq_aligner, gff] }
            .concat(BRAKER2.out.gff.map { id, meta, gff -> [id, meta, "braker2", gff] })
            .concat(GALBA.out.gff.map { id, meta, gff -> [id,  meta, "galba", gff] })
            .set { all_gff }

        // filtering annotation
        AGAT_STATISTICS(all_gff)
        AGAT_FILTER_GENE_BY_LENGTH(all_gff.map { id, meta, annotation_method, gff -> [id, meta, annotation_method, gff, 30] })

        AGAT_FILTER_GENE_BY_LENGTH.out
            .combine(masked_genomes, by: [0, 1])
            .map { id, meta, annotation_method, gff, genome -> [id, meta, annotation_method, genome, gff] }
            .set { gff_and_genome }

        EXTRACT_CANONICAL_PROTEINS(gff_and_genome)
        AGAT_FILTER_INCOMPLETE_GENES(gff_and_genome)

        DOWNLOAD_OMAMER_DB()
        EXTRACT_CANONICAL_PROTEINS.out.fasta
            .combine(DOWNLOAD_OMAMER_DB.out)
            .set { omark_input }

        // Run OMARK and BUSCO
        OMARK(omark_input)
        BUSCO(EXTRACT_CANONICAL_PROTEINS.out.fasta)

        OMARK.out.contaminants_count
            .filter { id, meta, annotation_method, contaminants_count -> contaminants_count.readLines()[0].toInteger() > 0 }
            .join ( EXTRACT_CANONICAL_PROTEINS.out.gff, by: [0, 1, 2] )
            .join ( EXTRACT_CANONICAL_PROTEINS.out.fasta, by: [0, 1, 2] )
            .join ( OMARK.out.ump, by: [0, 1, 2] )
            .map { id, meta, annotation_method, contaminants_count, gff, fasta, ump -> [id, meta, annotation_method, gff, fasta, ump, "0.33" ] }
            .set { contamination_filtering_input }

        contamination_filtering_input
            .concat (contamination_filtering_input.map{id, meta, annotation_method, gff, fasta, ump, threshold -> [id, meta, annotation_method, gff, fasta, ump, "0.1"]})
            .map {id, meta, annotation_method, gff, fasta, ump, threshold -> [ id, meta, annotation_method + "_filtered_c_"+ threshold, gff, fasta, ump, threshold ] }
            .set{contamination_filtering_input}

       // Check if OMARK filters enough contaminations
        FILTER_OMARK_CONTAMINATIONS(contamination_filtering_input)
        BUSCO_ON_CONTAMINATION_FILTERED(FILTER_OMARK_CONTAMINATIONS.out.fasta)
        OMARK_ON_CONTAMINATION_FILTERED(FILTER_OMARK_CONTAMINATIONS.out.fasta.map{id, meta, annotation_method, proteins -> [id, meta, annotation_method, proteins, luca_db] })

        if (params.module_repeats == true) {
            // Analysis overlapping between genes annotations and repeats
            EXTRACT_CANONICAL_PROTEINS.out.gff
                .join(repeats_gff, by: [0, 1]) // -> id, meta, annotation_method, proteins_gff, repeat_gff
                .set { genes_vs_repeats_input }

            CLASSIFY_GENES_VS_REPEATS(genes_vs_repeats_input)
        }
    // emit:
}


