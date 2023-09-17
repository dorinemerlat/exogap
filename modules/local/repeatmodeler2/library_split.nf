process LIBRARY_SPLIT {
    tag "LIBRARY_SPLIT_$library_id"
    cpus 32
    time '20d'

    input:
    tuple val(library_id), path(library_path)


    output:
    tuple val(library_id) path("${genome_id}-families-classified.fa"), emit: classified
    tuple val(library_id) path("${genome_id}-unclassified.fa"),        emit: unclassifier

    script:
    """
    sorting_repeats.py -i $library_path -k "${genome_id}-classified.fa" -u "${genome_id}-unclassified.fa"
    """
}
