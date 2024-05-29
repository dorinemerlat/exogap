process DOWNLOAD_NEWICK {
    // publishDir "${params.outdir}"
    label 'python_phylo'
    tag "newick"
    scratch false

    input:
    path(info)

    output:
    path "*.nwk"

    script:
    """
    get_newick.py $info
    """

    stub:
    """
    touch stub.nwk
    """
}
