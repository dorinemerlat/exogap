process AGAT_FILTER_INCOMPLETE_GENES {
    tag "${annotation_method}/${id}"
    // label 'agat'
    scratch false
    stageInMode 'copy'
    
    input:
    tuple val(id), val(meta), val(annotation_method), path(genome), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_complete.gff")

    script:
    """
    module load agat
    agat_sp_filter_incomplete_gene_coding_models.pl --gff $gff --fasta $genome -o ${id}_${annotation_method}_complete.gff
    """

    stub:
    """
    touch ${id}_${annotation_method}_complete.gff
    """
}