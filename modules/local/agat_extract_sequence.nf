process AGAT_EXTRACT_SEQUENCE {
    scratch true
    tag "AGAT_EXTRACT_SEQUENCE_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("proteins.gff"), val(iteration)

    script:
    """
    agat_sp_separate_by_record_type.pl -g $gff -o proteins.fa -f $genome --aa
    """

    stub:
    """
    touch proteins.gff
    """
}
