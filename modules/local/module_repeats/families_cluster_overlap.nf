process FAMILIES_CLUSTER_OVERLAP {
    tag "${id}"
    // label 'exogap_tools'

    input:
    tuple val(id), val(meta), path(cluster)

    output:
    tuple val(id), val(meta), path("${id}_parsed.tsv"), emit: tsv
    tuple val(id), val(meta), path("${id}_filtered.tsv"), emit: filtered_tsv
    tuple val(id), val(meta), path("${id}_summary_A4.png"), emit: summary_png

    script:
    """
    /home/merlat/.conda/envs/python/bin/python /tempor/merlat/exogaptwo/bin/families_cluster_overlap.py -n "$meta.name" -i $cluster -o $id
    """

    stub:
    """
    touch ${id}_parsed.tsv ${id}_filtered.tsv ${id}_summary_A4.png
    """
}
