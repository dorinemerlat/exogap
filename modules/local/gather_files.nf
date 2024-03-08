process gather_files {
    tag "gather_files_${name_output}"

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
}
