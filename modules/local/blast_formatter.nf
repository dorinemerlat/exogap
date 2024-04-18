process BLAST_FORMATTER {
    scratch true
    tag "BLAST_FORMATTER_${id}_${comment}"
    label 'blast'

    input:
    tuple val(id), val(meta), path(archive), val(format), val(comment)

    output:
    tuple val(id), val(meta), path("${id}.blastp_format_11.out"), val(comment)

    script:
    """
    blast_formatter -archive $archive -outfmt $format -out ${id}.blast.out

    """

    stub:
    """
    touch ${id}.blast.out
    """
}
