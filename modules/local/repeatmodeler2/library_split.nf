process LIBRARY_SPLIT {
    tag "LIBRARY_SPLIT_$library_id"

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
