process SUMMARIZE_REPEATS {
    tag "${id}"
    label 'exogap_python'
    memory { 5.GB * task.attempt }
    maxRetries 5

    input:
    tuple val(id), val(meta), path(gff)

    output:
    tuple val(id), val(meta), path("${id}_RE_summary.tsv")

    script:
    """
    summarize_repeats.py -i $gff -o ${id}_RE_summary.tsv -s ${meta.assembly_size}
    """

    stub:
    """
    touch ${id}_RE_summary.tsv
    """ 
}
