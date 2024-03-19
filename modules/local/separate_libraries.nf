process SEPARATE_LIBRARIES {
    tag "SEPARATE_LIBRARIES_$id"
    label 'bioawk'

    input:
    tuple val(id), val(meta), path(library)


    output:
    tuple val("${id}_classified"), val(meta), path("${id}_classified.fa"),     emit: classified
    tuple val("${id}_unclassified"), val(meta), path("${id}_unclassified.fa"), emit: unclassified

    script:
    """
    bioawk -c fastx '\$name !~ /#Unknown\$/ {print ">"\$name; print \$seq}' $library |fold -w 60  > ${id}_classified.fa
    bioawk -c fastx '\$name ~ /#Unknown\$/ {print ">"\$name; print \$seq}' $library |fold -w 60 > ${id}_unclassified.fa
    """

    stub:
    """
    touch ${id}_classified.fa ${id}_unclassified.fa
    """
}
