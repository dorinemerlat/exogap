process FILTER_MAKER_AB_INITIO_PREDICTIONS {
    scratch true
    tag "FILTER_MAKER_AB_INITIO_PREDICTIONS_${id}_${iteration}"

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("${id}_filtered.gff"), val(iteration)

    script:
    """
    grep -P 'snap|augustus' $gff > ${id}_filtered.gff
    """

    stub:
    """
    touch ${id}_filtered.gff
    """
}
