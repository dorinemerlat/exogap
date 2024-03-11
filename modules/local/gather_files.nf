process GATHER_FILES {
    // in tag, replace characters: [, ] by nothing
    tag "GATHER_FILES_${id}"

    input:
    tuple val(id), path(inputs), val(name_output), val(header)

    output:
    tuple  val(id), path(name_output)

    script:
    """
    if [ $header == no ]; then
        cat $inputs > ${name_output}

    elif [ $header == yes ]; then
        head -n 1 ${inputs[0]} > ${name_output}
        tail -n +2 $inputs >> ${name_output}
    fi
    """

    stub:
    """
    touch ${name_output}
    """
}
