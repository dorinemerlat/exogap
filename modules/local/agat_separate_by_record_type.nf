process AGAT_SEPARATE_BY_RECORD_TYPE {
    scratch true
    tag "AGAT_SEPARATE_BY_RECORD_TYPE_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("mrna.gff"), val(iteration)

    script:
    """
    agat_sp_separate_by_record_type.pl -g $gff -o mrna.gff
    """

    stub:
    """
    touch mrna.gff
    """
}
