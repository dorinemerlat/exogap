process SELECT_TRANSCRIPTS_FROM_SRA {
    tag "SELECT_TRANSCRIPTS_FROM_SRA_${id}"

    input:
    tuple val(id), val(meta), path(genome)

    output:
    stdout

    script:
    """

    """
}
