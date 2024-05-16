process DOWNLOAD_NEWICK {
    // publishDir "${params.outdir}"
    label 'python_phylo'
    tag "DOWNLOAD_NEWICK"

    input:
    path info

    output:
    path "*.tree"

    script:
    """
    get_newick.py $info
    """

    stub:
    """
    touch stub.tree
    """
}
