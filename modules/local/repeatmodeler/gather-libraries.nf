process GATHER_LIBRARIES {
    input:
    path library_paths

    output:
    tuple val("all-repeats"), path("all-repeats.fa")

    script:
    """
    cat $library_paths > all-repeats.fa
    """
}
