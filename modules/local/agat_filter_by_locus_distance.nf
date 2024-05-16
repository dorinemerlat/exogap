process AGAT_FILTER_BY_LOCUS_DISTANCE {
    tag "AGAT_FILTER_BY_LOCUS_DISTANCE_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("distance.gff"), val(iteration)

    script:
    """
    agat_sp_filter_by_locus_distance.pl --gff $gff -o distance.gff -d 1000
    """

    stub:
    """
    touch distance.gff
    """
}
