process LIBRARY_CONCATENATION {
    tag "LIBRARY_CONCATENATION_$genome_id"
    cpus 1

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.3--pl5321h9ee0642_0':
        'quay.io/biocontainers/repeatmodeler:2.0.3--pl5321h9ee0642_0' }"

    input:
    path library_paths

    output:
    tuple val("all_repeats_families"), path("all_repeats_families.fa")

    script:
    """
    cat $library_paths > all_repeats_families.fa
    """
}
