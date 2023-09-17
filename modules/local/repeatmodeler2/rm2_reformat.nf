process RM2_REFORMAT {
    tag "RM2_REFORMAT_$genome_id"

    input:
    tuple val(genome_id), path(genome_path)

    output:
    tuple val(genome_id), path("${genome_id}-families-renamed.fa")

    script:
    """
    rename_repeats.py.py -i g$enome_path -o ${genome_id}-families-renamed.fa -p repeats
    """
}
