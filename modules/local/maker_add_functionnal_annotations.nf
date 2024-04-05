process MAKER_ADD_FUNCTIONNAL_ANNOTATIONS {
    tag "MAKER_ADD_FUNCTIONNAL_ANNOTATIONS_${id}_${iteration}"
    label 'maker'
    cpus 40

    input:
    tuple val(id), val(meta), path(genome), path(train), path(test), path(specie), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path(train), path(test), path("config/species/${id}_*"), val(iteration)

    script:
    """
    cp -r /usr/local/config .
    export AUGUSTUS_CONFIG_PATH=\$(pwd)/config
    cp $specie config/species/

    optimize_augustus.pl --species=${id} --cpus=$task.cpus --kfold=$task.cpus --UTR=on $train
    """

    stub:
    """
    mkdir -p config/species/
    for i in $specie; do
        cp \$i config/species/
    done
    """
}
