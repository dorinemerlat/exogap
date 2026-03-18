process BUSCO {
    tag "all"
    label 'busco'
    scratch false
    cpus 10
    stageInMode 'copy'

    input:
    tuple val(id), val(json, stageAs: "input/*")

    output:
    tuple val(id), val(meta), path("busco_figure*png")

    script:
    """
    busco --plot input --plot_percentages
    """

    stub:
    """
    touch busco_figure.png
    """
}