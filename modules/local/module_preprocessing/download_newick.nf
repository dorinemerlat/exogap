process DOWNLOAD_NEWICK {
    tag "newick"
    label 'exogap_tools'

    input:
    tuple val(taxids)

    output:
    path "tree.nwk"

    script:
    """
    download_newick.py  --outfile tree.nwk $taxids
    """

    stub:
    """
    tree stub.nwk
    """
}
