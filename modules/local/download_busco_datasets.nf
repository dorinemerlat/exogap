process download_busco_datasets {
    tag "download_busco_datasets"
    cache 'lenient'
    label 'busco'

    output:
    path "BUSCO_datasets.txt"

    script:
    """
    busco --list-datasets | grep "_odb" |  awk '{print \$NF}' |awk -F '_' '{OFS=","; print \$1, \$1 "_" \$2}' > BUSCO_datasets.txt
    """
}
