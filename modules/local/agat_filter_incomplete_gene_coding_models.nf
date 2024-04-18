process AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS {
    scratch true
    tag "AGAT_FILTER_INCOMPLETE_GENE_CODING_MODELS_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("complete.gff"), val(iteration)

    script:
    """
    agat_sp_filter_incomplete_gene_coding_models.pl --gff $gff -o complete.gff -f $genome
    """

    stub:
    """
    touch complete.gff
    """
}
