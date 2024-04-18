process AGAT_FILTER_BY_MRNA_BLAST_VALUE {
    scratch true
    tag "AGAT_FILTER_BY_MRNA_BLAST_VALUE_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff), path(blastp_out), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("non_redundant.gff"), val(iteration)

    script:
    """
    agat_sp_filter_by_mrnaBlastValue.pl --gff $gff -o non_redundant.gff --blast $blastp_out
    """

    stub:
    """
    touch non_redundant.gff
    """
}
