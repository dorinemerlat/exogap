process DOWNLOAD_LINEAGE {
    tag "DOWNLOAD_LINEAGE_${id}"
    cache 'lenient'
    label 'jq'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), path("${id}.lineage"),   emit: lineage
    path "${id}_lineage.csv",               emit: lineage_for_info

    script:
    """
    # download the lineage and save it in a simple file to put it in the metadata
    curl -s https://lbgi.fr/api/taxonomy/lineage/$meta.taxid \
        | jq -r '.data[] | select(.rank) |"\\(.rank),\\(.id),\\(.name)"' \
        > ${id}.lineage

    # save it in another format to put it in the info file
    awk -F',' -v id=$id taxid=$meta.taxid 'BEGIN { FS=","; OFS="," } {
        keys[NR] = \$1
        values[NR] = \$2 "/" \$3
    }
    END {
        printf "id,taxid,"
        for (i = 1; i <= NR; i++) {
            printf "%s%s", keys[i], (i == NR ? "\\n" : ",")
        }
        printf id","taxid","
        for (i = 1; i <= NR; i++) {
            printf "%s%s", values[i], (i == NR ? "\\n" : ",")
        }
    }' ${id}.lineage > "${id}_lineage.csv
    """
}

