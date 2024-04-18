process DOWNLOAD_NEWICK {
    scratch true
    // publishDir "${params.outdir}"
    label 'python_phylo'
    tag "DOWNLOAD_NEWICK"

    input:
    tuple val(ids), val(meta), val(files)

    output:
    tuple val(ids), val(meta), path('*.tree'),             emit: newick
    tuple val(ids), val(meta), path('*.tree.ascii_art'),   emit: ascii

    script:
    """
    get-newick.py -i '$meta'
    """

    stub:
    """
    touch ${meta}.tree
    touch ${meta}.tree.ascii_art
    """
}
