process GATHER_FILES {
    scratch true
    // in tag, replace characters: [, ] by nothing
    tag "GATHER_FILES_${meta.name}"

    input:
    tuple val(id), val(meta), path(inputs), val(name_output), val(header)

    output:
    tuple  val(id), val(meta), path(name_output)

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
