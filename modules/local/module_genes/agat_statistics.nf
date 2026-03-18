process AGAT_STATISTICS {
    tag "${annotation_method}/${id}"
    stageInMode 'copy'
    // label 'agat'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}.stat")

    script:
    """
    module load agat
    agat_sp_statistics.pl --gff $gff  -o ${id}_${annotation_method}.stat
    """

    stub:
    """
    touch ${id}_${annotation_method}.stat
    """
}