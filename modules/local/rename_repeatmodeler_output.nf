process RENAME_REPEATMODELER_OUTPUT {
    tag "RENAME_REPEATMODELER_OUTPUT_${meta.id}"
    label 'biopython'

    input:
    tuple val(genome_id), path(genome_path)

    output:
    tuple val("${genome_id}"), path("${genome_id}-renamed.fa")

    script:
    """
    prefix=\$(echo $genome_id |sed -e "s/-families//g")
    rename_repeats.py -i $genome_path -o ${genome_id}-renamed.fa -p \$prefix
    """

    stub:
    """
    touch ${genome_id}-renamed.fa
    """
}

