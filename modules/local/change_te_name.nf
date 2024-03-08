process change_te_name {
    tag "change_te_name_${meta.id}"

    input:
    tuple val(meta), path(out)
    tuple val(meta), path(te_classes)

    output:
    tuple val(meta), path("${meta.id}_TE_found.csv"),   emit: publishing

    script:
    """
    awk -F ',' 'NR>1{print $0}' $te_classes | cut -f1,6,7 |awk '{OFS=","; print $1, $2 "/" $3}' |sed "s|/$||g" > ${te_class}.extract
    """
}
