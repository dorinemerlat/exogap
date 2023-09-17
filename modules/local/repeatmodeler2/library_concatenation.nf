process LIBRARY_CONCATENATION {
    input:
    path library_paths

    output:
    tuple val("all_repeats_families"), path("all_repeats_families.fa")

    script:
    """
    cat $library_paths > all_repeats_families.fa
    """
}
