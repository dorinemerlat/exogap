process GATHER_FAMILIES {
    // in tag, replace characters: [, ] by nothing
    tag "${id}"

    input:
    tuple val(id), val(meta), path(inputs)

    output:
    tuple val(id), val(meta), path("${id}_all_families.fa")

    script:
    """

    cat $inputs > ${id}_all_families.fa
    """

    stub:
    """
    touch ${id}_all_families.fa
    """
}
