process GENERATE_MAKER_GFF {
    tag "GENERATE_MAKER_GFF_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(maker_log), val(iteration)

    output:
    tuple  val(id), val(meta), path(genome), path("${id}_complete.gff"), val(iteration),                emit: gff_complete
    tuple  val(id), val(meta), path(genome), path("${id}.gff"), val(iteration),                         emit: gff
    tuple  val(id), val(meta), path(genome), path("${id}_complete_with_sequences.gff"), val(iteration), emit: gff_complete_with_sequences

    script:
    """
    gff3_merge -n -s -d $maker_log > ${id}_complete.gff # without sequences
    gff3_merge -s -d $maker_log > ${id}_complete_with_sequences.gff # with sequences

    grep -Pv '\tmatch\t|\tmatch_part\t|\tprotein_match\t|\texpressed_sequence_match\t' ${id}_complete.gff > ${id}.gff
    """

    stub:
    """
    touch ${id}.gff ${id}_complete_with_sequences.gff ${id}_complete.gff ${id}.gff
    """
}
