process AGAT_KEEP_LONGEST_ISOFORM {
    scratch true
    tag "AGAT_KEEP_LONGEST_ISOFORM_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("longest_isoform.gff"), val(iteration)

    script:
    """
    agat_sp_keep_longest_isoform.pl --gff $gff -o longest_isoform.gff
    """

    stub:
    """
    touch longest_isoform.gff
    """
}
