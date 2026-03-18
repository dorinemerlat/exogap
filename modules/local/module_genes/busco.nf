process BUSCO {
    tag "${annotation_method}/${id}"
    label 'busco'
    scratch false
    cpus 10
    stageInMode 'copy'

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins)

    output:
    tuple val(id), val(meta), val(annotation_method), path("short_summary.specific.arthropoda_odb10.busco_${id}_${annotation_method}.json")
    tuple val(id), val(meta), val(annotation_method), path("short_summary.specific.arthropoda_odb10.busco_${id}_${annotation_method}.txt")

    script:
    """
    busco -i $proteins -o busco_${id}_${annotation_method} -m prot -l arthropoda_odb10 -c $task.cpus -f
    mv busco_${id}_${annotation_method}/short_summary*.{json,txt} .
    """

    stub:
    """
    touch short_summary.specific.arthropoda_odb10.busco_${id}_${annotation_method}.json
    touch short_summary.specific.arthropoda_odb10.busco_${id}_${annotation_method}.txt
    """
}