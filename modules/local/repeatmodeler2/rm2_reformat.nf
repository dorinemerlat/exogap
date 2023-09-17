process RM2_REFORMAT {
    tag "RM2_REFORMAT_$genome_id"

    input:
    tuple val(genome_id), path(genome_path)

    output:
    tuple val("${genome_id}-renamed"), path("${genome_id}-renamed.fa")

    script:
    """
    prefix=\$(echo $genome_id |sed -e "s/-families//g")
    rename_repeats.py -i $genome_path -o ${genome_id}-renamed.fa -p \$prefix
    """
}
