process RENAME_GENOME {
    tag "RENAME_GENOME_$id"
    // publishDir "${params.outdir}/out/$id"
    cache 'lenient'
    label 'bioawk'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path(genome),              emit: fasta
    tuple val(id), val(meta), path("${id}.tsv"),         emit: json

    script:
    def unmodified_genome = "${id}_unreformated.fa"
    """
    # Name modification to avoid overwriting the original file (symbolic link)
    mv $genome $unmodified_genome

    # Saving old names (with the new names) in a csv file
    bioawk  -v id=$id -c fastx 'BEGIN{print "previous_IDs", "EXOGAP_IDs"} {print \$name, id"_seq"NR}' $unmodified_genome > ${id}.tsv

    # Modification of genome file
    bioawk -v id=$id -c fastx '{print ">"id"_seq"NR; print \$seq}' $unmodified_genome |fold -w 60 > $genome
    """

    stub:
    """
    touch $genome
    touch ${id}.tsv
    """
}
