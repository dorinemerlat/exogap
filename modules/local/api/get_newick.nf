process GET_NEWICK {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://dorinemerlat/python-exogap:v1.02':
    'dorinemerlat/python-exogap:v1.02' }"

    cache 'lenient'

    publishDir "${params.outdir}/out/"

    input:
    tuple val(ids), val(metas), val(files)

    output:
    tuple val(ids), val(metas), path('*.tree'),             emit: newick
    tuple val(ids), val(metas), path('*.tree.ascii_art'),   emit: ascii

    script:
    """
    get-newick.py -i '$metas'
    """

    stub:
    """
    touch ${metas}.tree
    touch ${metas}.tree.ascii_art
    """
}
