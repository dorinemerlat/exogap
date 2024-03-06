process SELECT_PUBLIC_DATA {
    tag "SELECT_PUBLIC_DATA_${id}"
    publishDir "${params.outdir}/results/data_to_use"

    input:
    tuple val(id), val(meta), path(genome), val(taxons), path(info)

    output:
    tuple val(id), val(meta), path(genome), path("${id}_data_used.infos"), path("${id}_data_used.infos.temp")

    script:
    def max_proteins_from_proteomes = params.max_proteins_from_a_large_set_of_species
    def max_proteins_from_relative  = params.max_proteins_from_relatively_close_species
    def max_proteins_from_close     = params.max_proteins_from_close_species
    def max_transcriptomes          = params.max_transcriptomes


    """
    for taxon in \$(echo ${taxons}); do
        echo ${id},\${taxon} >> ${id}_infos.csv
    done

    join -t "," -1 2 -2 1 -o 1.1,1.2,2.2,2.3,2.4,2.5 <(sort -t, -k2 ${id}_infos.csv) <(sort -t, -k1 ${info}) > infos.tsv

    awk -F, -v max=${max_proteins_from_relative} '\$4 < max' infos.tsv | sort -t ',' -k4,4nr |head -n 1 > infos_max_proteins_from_relative.tsv
    awk -F, -v max=${max_proteins_from_close} '\$4 < max' infos.tsv | sort -t ',' -k4,4nr |head -n 1 > infos_max_proteins_from_close.tsv
    awk -F, -v max=${max_proteins_from_proteomes} '\$5 < max' infos.tsv | sort -t ',' -k5,5nr |head -n 1 > infos_max_proteins_from_proteomes.tsv
    awk -F, -v max=${max_transcriptomes} '\$6 < max' infos.tsv | sort -t ',' -k6,6nr |head -n 1 > infos_max_transcriptomes.tsv

    echo "Data used for annotation of ${id}: \n" > ${id}_data_used.infos
    awk -F, '{print "Proteins from a set of large set of species:       "\$5" protein sequences from "\$3" ("\$2")"}' infos_max_proteins_from_proteomes.tsv >> ${id}_data_used.infos
    awk -F, '{print "Proteins from a set of very close species:         "\$4" protein sequences from "\$3" ("\$2")"}' infos_max_proteins_from_close.tsv >> ${id}_data_used.infos
    awk -F, '{print "Proteins from a set of relatively close species:   "\$4" protein sequences from "\$3" ("\$2")"}' infos_max_proteins_from_relative.tsv >> ${id}_data_used.infos
    awk -F, '{print "Transcripts from set of species:                   "\$6" transcriptomes from "\$3" ("\$2")"}' infos_max_transcriptomes.tsv >> ${id}_data_used.infos

    awk -F, '{print "large_protein_set,"\$2}' infos_max_proteins_from_proteomes.tsv > ${id}_data_used.infos.temp
    awk -F, '{print "relatively_close_protein_set,"\$2}' infos_max_proteins_from_relative.tsv >> ${id}_data_used.infos.temp
    awk -F, '{print "very_close_protein_set,"\$2}' infos_max_proteins_from_close.tsv >> ${id}_data_used.infos.temp
    awk -F, '{print "transcriptome_set,"\$2}' infos_max_transcriptomes.tsv >> ${id}_data_used.infos.temp
    """
}

    // echo meta.lineage

    // genome_id,taxon_taxid

    // -->
    // genome_id,taxon_taxid,taxon_name,proteins_count,proteins_from_proteomes_count,transcriptomes_count
