process GENERATE_MAKER_FASTA {
    tag "GENERATE_MAKER_FASTA_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(maker_log), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("${id}.proteins.fasta"),       emit: proteins
    tuple val(id), val(meta), path(genome), path("${id}.transcripts.fasta"),    emit: transcripts


    script:
    """
    fasta_merge -d ${maker_log} -o ${id}
    sed -i 's/\\*//g' ${id}*.fasta
    """

    stub:
    """
    touch ${id}.proteins.fasta ${id}.transcripts.fasta
    """
}
