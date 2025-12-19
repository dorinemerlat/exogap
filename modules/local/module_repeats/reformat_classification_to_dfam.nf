process REFORMAT_CLASSIFICATION_TO_DFAM {
    tag "${id}"
    cpus 10
    memory { (10.GB * (task.attempt * task.attempt)) }
    maxRetries 5
    label 'exogap_tools'

    input:
    tuple val(id), val(meta), path(gff)
    path dfam

    output:
    tuple val(id), val(meta), path("${id}_repeats_reformated.gff")

    script:
    """
    reformat_classification_to_dfam.sh -g "$gff" -d "$dfam" -o "${id}_repeats_reformated.gff" -t $task.cpus
    """

    stub:
    """
    touch ${id}_repeats_reformated.gff
    """
}
