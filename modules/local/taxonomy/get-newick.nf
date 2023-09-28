process GET_NEWICK {
    debug
    publishDir "results/"
    // conda (params.enable_conda ? 'bioconda requests==2.26.0' : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/requests:2.26.0':
    //     'quay.io/biocontainers/requests:2.26.0' }"

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
