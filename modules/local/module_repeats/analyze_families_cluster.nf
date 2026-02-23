process ANALYZE_FAMILIES_CLUSTER {
    tag "${id}"
    // label 'exogap_tools'
    scratch 'false'

    input:
    tuple val(id), val(meta), path(cluster), val(out_prefix)

    output:
    tuple val(id), path("{out_prefix}_parsed.tsv"), emit: tsv
    tuple val(id), path("${out_prefix}_filtered.tsv"), emit: filtered_tsv
    tuple val(id), path("${out_prefix}_upset_TE_clusters.png"), emit: upset_plot
    tuple val(id), path("${out_prefix}_venn_TE_clusters.png"), emit: venn
    tuple val(id), path("${out_prefix}_cluster_size_distribution.png"), emit: cluster_size_dist
    tuple val(id), path("${out_prefix}_confirmation_pies.png"), emit: confirmation_p

    script:
    """
    scluster_species_overlap.py -i $cluster -o $out_prefix
    """

    stub:
    """
    touch ${out_prefix}_parsed.tsv ${out_prefix}_filtered.tsv ${out_prefix}_upset_TE_clusters.png ${out_prefix}_venn_TE_clusters.png \
        ${out_prefix}_cluster_size_distribution.png ${out_prefix}_confirmation_pies.png
    """
}
