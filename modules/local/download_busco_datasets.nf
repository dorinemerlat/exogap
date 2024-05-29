process DOWNLOAD_BUSCO_DATASETS {
    tag "datasets"
    cache 'lenient'
    label 'busco'
    scratch false 

    input:
    tuple val(id), val(list_datasets)

    output:
    path "busco_downloads"

    script:
    """
    busco --download $list_datasets
    busco --download \$(awk '\$NF == "placement_files" {print \$1}' busco_downloads/file_versions.tsv | tr '\n' ' ')
    """
}
