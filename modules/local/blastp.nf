process BLASTP {
    tag "BLASTP_${id}_${comment}"
    label 'blast'

    input:
    tuple val(id), val(meta), path(proteins), val(db), path(db_files), val(evalue), val(word_size), val(max_hsps), val(format), val(comment)

    output:
    tuple val(id), val(meta), path("${id}.blastp.out"), val(comment)

    script:
    """
    blastp -query $query -db $db -evalue $evalue -word_size $word_size -max_hsps $max_hsps -outfmt $format -out ${id}.blastp.out
    """

    stub:
    """
    touch ${id}.blastp.out
    """
}
