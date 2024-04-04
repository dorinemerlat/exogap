process BLASTP {
    tag "BLASTP_${id}_${comment}"
    label 'blast'

    input:
    tuple val(id), val(meta), path(proteins), path(db), val(format), val(comment)

    output:
    tuple val(id), val(meta), path("${id}.blastp.out"), val(comment)

    script:
    """
    blastp -query $query -db $db -outfmt $format -out ${id}.blastp.out
    """

    stub:
    """
    touch ${id}.blastp.out
    """
}
