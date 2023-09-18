process UNCLASSIFIED_SEPARATION {
    tag "UNCLASSIFIED_SEPARATION_$library_id"

    conda (params.enable_conda ? 'bioawk==1.0--h7132678_8' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h7132678_7':
        'quay.io/biocontainers/bioawk:1.0--h7132678_8' }"

    input:
    tuple val(library_id), path(library_path)


    output:
    tuple val("${library_id}-classified"), path("${library_id}-classified.fa"), emit: classified
    tuple val("${library_id}-unclassified"), path("${library_id}-unclassified.fa"),        emit: unclassified

    script:
    """
    sorting_repeats.py -i $library_path -c "${library_id}-classified.fa" -u "${library_id}-unclassified.fa"
    """
}
