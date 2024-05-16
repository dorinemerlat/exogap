process MAKER_MAP_IDS {
    tag "MAKER_MAP_IDS_${id}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(gff), path(proteins), path(transcripts)

    output:
    tuple val(id), val(meta), path("${gff}.id.map"),                  emit: map
    tuple val(id), val(meta), path("${id}_renamed.gff"),              emit: gff
    tuple val(id), val(meta), path("${id}_renamed.proteins.fa"),      emit: proteins
    tuple val(id), val(meta), path("${id}_renamed.transcripts.fa"),   emit: transcripts

    script:
    """
    cp $gff ${id}_renamed.gff
    cp $proteins ${id}_renamed.proteins.fa
    cp $transcripts ${id}_renamed.transcripts.fa

    maker_map_ids --prefix ${id}_ $gff > ${gff}.id.map
    map_gff_ids ${gff}.id.map ${id}_renamed.gff
    map_fasta_ids  ${gff}.id.map ${id}_renamed.proteins.fa
    map_fasta_ids  ${gff}.id.map ${id}_renamed.transcripts.fa
    """

    stub:
    """
    touch ${gff}.id.map ${id}_renamed.gff ${id}_renamed.proteins.fa ${id}_renamed.transcripts.fa
    """
}
