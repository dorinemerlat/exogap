process GET_TAXONOMIC_LINEAGE {
    tag "GET_TAXONOMIC_LINEAGE_${id}"

    input:
    tuple val(id), val(meta), path(genome)

    output:
    stdout

    script:
    """
    curl -s https://lbgi.fr/api/taxonomy/lineage/$meta.taxid | jq -r '.data[] | select(.rank) |"\\(.rank),\\(.id),\\(.name)"' \
        | awk -v id=$id 'BEGIN{print id} {print \$0}'
    """
}

