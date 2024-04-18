process GFF_TO_AUGUSTUS_GFF {
    scratch true
    tag "GFF_TO_AUGUSTUS_GFF_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple  val(id), val(meta), path(genome), path("${id}_for_augustus.gff"), val(iteration)

    script:
    """
    python GFF2augustusGFF.py $gff ${id}_for_augustus.gff
    """

    stub:
    """
    touch ${id}_for_augustus.gff
    """
}
