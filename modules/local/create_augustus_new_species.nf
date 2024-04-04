process CREATE_AUGUSTUS_NEW_SPECIES {
    tag "CREATE_AUGUSTUS_NEW_SPECIES_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("${id}_*"), val(iteration)

    script:
    """
    cp -r /usr/local/config .
    export AUGUSTUS_CONFIG_PATH=\$(pwd)/config

    new_species.pl --species=${augustus_specie}
    """

    stub:
    """
    for i in ["metapars.cfg", "metapars.cgp.cfg", "metapars.utr.cfg", "parameters.cfg", "weightmatrix.txt"]; do
        touch ${id}_\${i}
    done
    """
}
