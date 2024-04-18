process DOWNLOAD_LINEAGE {
    scratch true
    tag "DOWNLOAD_LINEAGE_${id}"
    cache 'lenient'
    label 'jq'

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

