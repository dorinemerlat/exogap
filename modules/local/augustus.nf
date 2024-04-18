process AUGUSTUS {
    scratch true
    tag "AUGUSTUS_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(train), path(test), path(specie), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("augustus_${iteration}.out"), val(iteration)

    script:
    """
    cp -r /usr/local/config .
    export AUGUSTUS_CONFIG_PATH=\$(pwd)/config
    cp $specie config/species/

    augustus --species=${id} $test | tee augustus_${iteration}.out
    """

    stub:
    """
    touch augustus_${iteration}.out
    """
}
