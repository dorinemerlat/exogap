process fasta_validator {
    tag "fasta_validator_${id}"
    cache 'lenient'
    label 'fasta_validator'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path(genome)

    script:
    """
    py_fasta_validator -f $genome
    """

    stub:
    """
    touch run_fasta_validator.stub
    """
}
