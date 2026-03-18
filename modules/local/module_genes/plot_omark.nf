process OMARK {
    tag "${annotation_method}/${id}"
    // label 'agat'
    scratch false
    memory '20 GB'
    stageInMode 'copy'
    
    input:
    tuple val(id), val(omark, stageAs: "input/*")

    output:
    tuple val(id), val(meta), path("omark"), val(annotation_method), emit: omark_out
    tuple val(id), val(meta), path("${id}_${annotation_method}_contaminations.txt"), val(annotation_method), emit: contaminations_txt
    tuple val(id), val(meta), path("${id}_${annotation_method}_contaminations.fa"), val(annotation_method), emit: contaminations_fa

    script:
    """
    plot_all_results.py -i input -o omark 
    """

    stub:
    """
    touch omark/${id}_${annotation_method}_detailed_summary.txt ${id}_${annotation_method}_contaminations.txt ${id}_${annotation_method}_contaminations.fa
    """
}
