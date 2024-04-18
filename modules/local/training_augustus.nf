process TRAINING_AUGUSTUS {
    scratch true
    tag "TRAINING_AUGUSTUS_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(train), path(test), path(specie), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path(train), path(test), path("config/species/${id}_*"), val(iteration)

    script:
    """
    cp -r /usr/local/config .
    export AUGUSTUS_CONFIG_PATH=\$(pwd)/config
    cp $specie config/species/

    etraining --species=${id} $train
    """

    stub:
    """
    mkdir -p config/species/
    for i in $specie; do
        cp \$i config/species/
    done
    """
}
