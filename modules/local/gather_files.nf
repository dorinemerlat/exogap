process GATHER_FILES {
    scratch true
    // in tag, replace characters: [, ] by nothing
    tag "GATHER_FILES_${meta.name}"

    input:
    tuple val(id), val(meta), path(inputs), val(output_name), val(header)

    output:
    tuple val(id), val(meta), path(output_name)

    script:
    """
    if [ $header == no ]; then
        cat $inputs > ${output_name}

    elif [ $header == yes ]; then
        head -n 1 ${inputs[0]} > ${output_name}
        tail -n +2 $inputs >> ${output_name}
    fi
    """

    stub:
    """
    touch ${output_name}
    """
}
