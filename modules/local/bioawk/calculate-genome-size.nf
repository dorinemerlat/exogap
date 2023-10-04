process CALCULATE_GENOME_SIZE {
    tag "CALCULATE_GENOME_SIZE_${meta.id}"

    conda (params.enable_conda ? 'bioawk==1.0--h7132678_8' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h7132678_7':
        'quay.io/biocontainers/bioawk:1.0--h7132678_8' }"

    input:
    tuple val(meta), path(genome)


    output:
    tuple val(meta), path("${meta.id}.csv")

    script:
    """
    bioawk -c fastx '{print length(\$seq)}' $genome |awk -v id=$meta.id 'BEGIN{OFS=",";} {sum+=\$1} END {print id, sum;}' > "${meta.id}.csv"
    """
}
