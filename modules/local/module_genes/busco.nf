process BUSCO {
    tag "${annotation_method}/${id}"
    label 'busco'
    scratch false
    cpus 10
    stageInMode 'copy'

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_busco.json")
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_busco.txt")

    script:
    """
    busco -i $proteins -o busco_${id}_${annotation_method} -m prot -l arthropoda_odb10 -c ${task.cpus} -f
    mv busco_${id}_${annotation_method}/short_summary*.json,txt busco_${id}_${annotation_method}_busco.json
    mv busco_${id}_${annotation_method}/short_summary*.txt busco_${id}_${annotation_method}_busco.txt
    """

    stub:
    """
    touch busco_${id}_${annotation_method}_busco.json
    touch busco_${id}_${annotation_method}_busco.txt
    """
}