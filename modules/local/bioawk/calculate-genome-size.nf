process CALCULATE_GENOME_SIZE {
    tag "CALCULATE_GENOME_SIZE_$library_id"

    conda (params.enable_conda ? 'bioawk==1.0--h7132678_8' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h7132678_7':
        'quay.io/biocontainers/bioawk:1.0--h7132678_8' }"

    input:
    tuple val(genome_id), path(genome_path)


    output:
    tuple val("${genome_id}"), path("${genome_id}-genome-size.csv")

    script:
    """
    bioawk -c fastx '{print length(\$seq)}' $genome_path |awk -v id=$genome_id 'BEGIN{OFS=",";} {sum+=\$1} END {print id, sum;}' > ${genome_id}-genome-size.csv
    """
}
