process RENAME_GENOMES {
    tag "RENAME_GENOMES_$id"
    publishDir "${params.outdir}/results/$id"

    conda (params.enable_conda ? 'bioawk==1.0--h7132678_8' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h7132678_7':
        'quay.io/biocontainers/bioawk:1.0--h7132678_8' }"

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path(genome),              emit: fasta
    tuple val(id), val(meta), path("${id}.tsv"),    emit: json

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
}
