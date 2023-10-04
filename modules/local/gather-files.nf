process GATHER_FILES {
    tag "GATHER_FILES_${name}"

    input:
    tuple val(name), path(all_files)

    output:
    tuple val(name), path("all-${name}")

    script:
    """
    cat $all_files > "all-${name}"
    """
}
