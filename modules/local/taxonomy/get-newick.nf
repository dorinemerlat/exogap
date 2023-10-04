process GET_NEWICK {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://dorinemerlat/python-exogap:v1.01':
    'dorinemerlat/python-exogap:v1.01' }"

    publishDir "${params.outdir}/results/"

    input:
    val(IDs)

    output:
    path('*.tree'),             emit: newick
    path('*.tree.ascii_art'),   emit: ascii

    script:
    """
    get-newick.py -i '$IDs'
    """
}
