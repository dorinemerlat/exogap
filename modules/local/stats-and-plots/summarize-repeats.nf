process SUMMARIZE_REPEATS {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://dorinemerlat/python-exogap:v1.02':
    'dorinemerlat/python-exogap:v1.02' }"

    input:
    path(out)
    path(classification)
    path(newick)

    output:
    path("TE-overview.pdf")

    script:
    """
    touch TE-overview.pdf
    # summarize-repeats.py -c $classification -n $newick
    """
}
