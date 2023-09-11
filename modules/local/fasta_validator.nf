process FASTA_VALIDATOR {
    tag "FASTA_VALIDATOR_$genome_id"

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/py_fasta_validator:0.5--py39hc16433a_3':
        'quay.io/biocontainers/py_fasta_validator:0.5--py39h7d875b9_0' }"

    input:
    tuple val(genome_id), path(genome_path)

    output:
    tuple val(genome_id), path(genome_path)

    script:
    """
    py_fasta_validator -f $genome_path
    """
}
