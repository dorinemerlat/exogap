process EXTRACT_SEQUENCES {
    scratch true
    tag "EXTRACT_SEQUENCES_${meta.id}"
    publishDir "${params.outdir}/results/${meta.id}/repetitive-elements"

    label 'agat'

    input:
    tuple val(meta), path(genome), val(type)
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${genome}-repetitive-elements.fa")

    script:
    """
    agat_sp_extract_sequences.pl -g $gff -f $genome -t $type -o ${genome}-repetitive-elements.fa
    """
}

// agat_sp_extract_sequences.pl -g $gff -f $genome -t dispersed_repeat -o ${genome}-repetitive-elements.fa
