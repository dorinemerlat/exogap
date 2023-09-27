process GATHER_FILES {
    debug

    input:
    tuple val(file_id), path(file_path)

    output:
    tuple val(file_id), path("${file_id}.all")

    script:
    """
    cat $file_path > ${file_id}.all
    """
}
