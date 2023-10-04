process FASTA_VALIDATOR {
    tag "FASTA_VALIDATOR_${meta.id}"

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/py_fasta_validator:0.5--py39hc16433a_3':
        'quay.io/biocontainers/py_fasta_validator:0.5--py39h7d875b9_0' }"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path(genome)

    script:
    """
    py_fasta_validator -f $genome
    """
}
