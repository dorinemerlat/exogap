process bioawk_to_separate_repeat_libraries {
    tag "bioawk_to_separate_repeat_libraries_$library_id"
    label 'bioawk'

    input:
    tuple val(library_id), path(library_path)


    output:
    tuple val("${library_id}-classified"), path("${library_id}-classified.fa"),            emit: classified
    tuple val("${library_id}-unclassified"), path("${library_id}-unclassified.fa"),        emit: unclassified

    script:
    """
    bioawk -c fastx '\$name !~ /#Unknown\$/ {print ">"\$name; print \$seq}' $library_path |fold -w 60  > ${library_id}-classified.fa
    bioawk -c fastx '\$name ~ /#Unknown\$/ {print ">"\$name; print \$seq}' $library_path |fold -w 60 > ${library_id}-unclassified.fa
    """
}
