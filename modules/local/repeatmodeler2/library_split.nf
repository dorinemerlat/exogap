process LIBRARY_SPLIT {
    tag "LIBRARY_SPLIT_$genome_id"
    cpus 32
    time '20d'

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.3--pl5321h9ee0642_0':
        'quay.io/biocontainers/repeatmodeler:2.0.3--pl5321h9ee0642_0' }"

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
