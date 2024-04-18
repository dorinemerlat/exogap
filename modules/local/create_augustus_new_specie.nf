process CREATE_AUGUSTUS_NEW_SPECIE {
    scratch true
    tag "CREATE_AUGUSTUS_NEW_SPECIE_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("config/species/${id}_*"), val(iteration)

    script:
    """
    cp -r /usr/local/config .
    export AUGUSTUS_CONFIG_PATH=\$(pwd)/config

    new_species.pl --species=${id}
    """

    stub:
    """
    mkdir -p config/species
    for i in "metapars.cfg", "metapars.cgp.cfg", "metapars.utr.cfg", "parameters.cfg", "weightmatrix.txt"; do
        touch config/species/${id}_\${i}
    done
    """
}
