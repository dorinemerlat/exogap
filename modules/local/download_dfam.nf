process download_dfam {
    output:
    path("dfam_classification.csv")

    script:
    """
    curl -s https://dfam.org/api/classes \
        | jq -r '. as $parent | .. | select(.full_name?) | [.full_name, .repeatmasker_type, .repeatmasker_subtype] | @tsv' \
        | sed 's/^root;//g' \
        | awk '{OFS=","; print $1, $2 "/" $3}' \
        | sed "s|/$||g" \
        > dfam_classification.csv
    """
}
