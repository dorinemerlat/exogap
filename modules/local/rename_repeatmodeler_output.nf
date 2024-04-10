process RENAME_REPEATMODELER_OUTPUT {
    tag "RENAME_REPEATMODELER_OUTPUT_${id}"
    label 'biopython'

    input:
    tuple  val(id), val(meta), path(repeatmodeler_output)

    output:
    tuple  val(id), val(meta), path("${id}_renamed.fa")

    script:
    """
    prefix=\$(echo $repeatmodeler_output |sed -e "s/-families//g")
    rename_repeats.py -i $repeatmodeler_output -o ${genome_id}_renamed.fa -p \$prefix
    """

    stub:
    """
    touch ${id}_renamed.fa
    """
}

