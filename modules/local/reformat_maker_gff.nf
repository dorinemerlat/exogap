process REFORMAT_MAKER_GFF {
    tag "REFORMAT_MAKER_GFF_${id}_${iteration}"
    label 'biopython'

    input:
    tuple val(id), val(meta), path(genome), path(gff) // query

    output:
    tuple val(id), val(meta), path(genome), path("${id}_final.gff") // query
    script:
    """
    touch ${id}_final.gff
    """

    stub:
    """
    touch ${id}_final.gff
    """
}
