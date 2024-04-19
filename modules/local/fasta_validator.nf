process FASTA_VALIDATOR {
    tag "FASTA_VALIDATOR_${id}"
    cache 'lenient'
    label 'fasta_validator'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path(genome), path('run_fasta_validator.*'),   emit: fasta

    script:
    """
    py_fasta_validator -f $genome
    echo '${genome} is a valid fasta' > run_fasta_validator.out
    """

    stub:
    """
    touch run_fasta_validator.out
    """
}
