process GATHER_FILES {
    tag "GATHER_FILES_${name_output}"

    input:
    tuple val(id), path(inputs), val(name_output), val(ext_output), val(header)

    output:
    tuple  val(id), path("all_${name_output}.${ext_output}")

    script:
    """
    if [ $header == no ]; then
        cat $inputs > all_${name_output}.${ext_output}

    elif [ $header == yes ]; then
        head -n 1 ${inputs[0]} > all_${name_output}.${ext_output}
        tail -n +2 $inputs >> all_${name_output}.${ext_output}
    fi
    """
}
