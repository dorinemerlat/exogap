process GENERATE_MAKER_GFF {
    tag "GENERATE_MAKER_GFF_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(maker_log), val(iteration)

    output:
    tuple  val(id), val(meta), path("${id}.gff"), val(iteration)
    tuple  val(id), val(meta), path("${id}_complete_with_sequences.gff"), val(iteration)

    script:
    """
    gff3_merge -n -s -d $maker_log > ${id}_complete.gff # without sequences
    gff3_merge -s -d $maker_log > ${id}_complete_with_sequences.gff # with sequences

    grep -Pv '\tmatch\t|\tmatch_part\t|\tprotein_match\t|\texpressed_sequence_match\t' ${id}_complete.gff > ${id}.gff
    """

    stub:
    """
    touch ${id}.gff
    """
}
