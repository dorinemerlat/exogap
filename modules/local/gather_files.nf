process GATHER_FILES {
    // in tag, replace characters: [, ] by nothing
    tag "${id}"

    input:
    tuple val(id), val(meta), path(inputs), val(output_name), val(header)

    output:
    tuple val(id), val(meta), path(output_name)

    script:
    """
    if [ $header == no ]; then
        cat $inputs > ${output_name}

    elif [ $header == yes ]; then
        first_input=\$(echo "$inputs" | cut -f 1 -d ' ')
        head -n 1 \$first_input > ${output_name}
        tail -q -n +2 $inputs >> ${output_name}
    fi
    """

    stub:
    """
    touch ${output_name}
    """
}