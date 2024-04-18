process RENAME_GENOME {
    scratch true
    tag "RENAME_GENOME_$id"
    // publishDir "${params.outdir}/out/$id"
    cache 'lenient'
    label 'bioawk'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("genome_${id}.fa"),      emit: fasta
    tuple val(id), val(meta), path("genome_${id}.fa.tsv"),  emit: json

    script:
    def unmodified_genome = "${id}_unreformated.fa"
    """
    # Name modification to avoid overwriting the original file (symbolic link)
    mv $genome $unmodified_genome

    # Saving old names (with the new names) in a csv file
    bioawk  -v id=$id -c fastx 'BEGIN{print "previous_IDs", "EXOGAP_IDs"} {print \$name, id"_seq"NR}' $unmodified_genome > genome_${id}.fa.tsv

    # Modification of genome file
    bioawk -v id=$id -c fastx '{print ">"id"_seq"NR; print \$seq}' $unmodified_genome |fold -w 60 > genome_${id}.fa
    """

    stub:
    """
    touch genome_${id}.fa
    touch genome_${id}.fa.tsv
    """
}
