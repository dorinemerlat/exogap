process download_lineage {
    tag "download_lineage_${id}"
    cache 'lenient'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), path("${id}.lineage")

    script:
    """
    curl -s https://lbgi.fr/api/taxonomy/lineage/$meta.taxid \
        | jq -r '.data[] | select(.rank) |"\\(.rank),\\(.id),\\(.name)"' \
        > ${id}.lineage
    """
}

