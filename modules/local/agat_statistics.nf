process AGAT_STATISTICS {
    tag "AGAT_STATISTICS_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration) // without sequences

    output:
    tuple val(id), val(meta), path(genome), path("${id}.stats"), val(iteration)

    script:
    """
    agat_sp_statistics.pl --gff $gff -o ${id}.stats
    """

    stub:
    """
    touch ${id}.stats
    """
}
