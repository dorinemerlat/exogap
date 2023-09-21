process GATHER_LIBRARIES {
    input:
    path library_paths

    output:
    tuple val("all-repeats-families"), path("all-repeats-families.fa")

    script:
    """
    more $library_paths > all-repeats-families.fa
    """
}
