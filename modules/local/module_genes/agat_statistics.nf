process AGAT_STATISTICS {
    tag "${annotation_method}/${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_stats.out")

    script:
    """
    agat_sp_statistics.pl --gff $gff  -o ${id}_${annotation_method}_stats.out
    """

    stub:
    """
    touch ${id}_${annotation_method}_stats.out
    """
}