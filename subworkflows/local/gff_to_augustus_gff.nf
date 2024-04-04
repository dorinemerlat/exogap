process GFF_TO_AUGUSTUS_GFF {
    tag "GFF_TO_AUGUSTUS_GFF_${id}_${iteration}"
    label 'biopython'

    input:
    tuple val(id), val(meta), path(genome), path(gff)
    each val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("augustus.gff")

    script:
    """
    GFF2augustusGFF.py gff augustus.gff
    """

    stub:
    """
    touch augustus.gff
    """
}
