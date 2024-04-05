process MAKEBLASTDB {
    tag "MAKEBLASTDB_${id}_${comment}"
    label 'blast'

    input:
    tuple val(id), val(meta), path(sequences), val(dbtype), val(comment)

    output:
    tuple val(id), val(meta), path("${id}_db.*"), val(comment)

    script:
    """
    makeblastdb -in $sequences -dbtype $dbtype -out ${id}_db
    """

    stub:
    """
    for i in "ndb", "nhr", "nin", "not", "nsq", "ntf", "nto"; do
        touch ${id}_db.\${i}
    done
    """
}
