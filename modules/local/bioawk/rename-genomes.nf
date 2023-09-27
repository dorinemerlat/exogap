process RENAME_GENOMES {
    tag "RENAME_GENOMES_$genome_id"
    publishDir "results/$genome_id"
    debug

    conda (params.enable_conda ? 'bioawk==1.0--h7132678_8' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h7132678_7':
        'quay.io/biocontainers/bioawk:1.0--h7132678_8' }"

    input:
    tuple val(genome_id), path(genome_path)

    output:
    tuple val(genome_id), path(genome_path),            emit: fasta
    tuple val(genome_id), path("${genome_id}.tsv"),    emit: json

    script:
    """
    mv $genome_path ${genome_id}_unreformated.fa
    grep '>' ${genome_id}_unreformated.fa | sed -e "s/^>//g" | awk -v id=$genome_id '{print \$O, id"-seq"NR}' > ${genome_id}.tsv
    bioawk -v id=$genome_id -c fastx '{print ">"id"-seq"NR; print \$seq}' ${genome_id}_unreformated.fa |fold -w 60 > $genome_path
    """
}
