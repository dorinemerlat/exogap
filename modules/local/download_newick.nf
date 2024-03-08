process download_newick {
    cache 'lenient'
    publishDir "${params.outdir}/out/"
    label 'python_phylo'

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
