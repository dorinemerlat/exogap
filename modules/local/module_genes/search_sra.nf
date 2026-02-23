process SEARCH_SRA {
    tag "${name}"
    label 'exogap_tools'
    // label 'retry_with_backoff'

    input:
    tuple val(dataset_taxid), val(dataset_name), val(taxid), val(name)

    output:
    tuple val(dataset_taxid), val(dataset_name), val(taxid), val(name), path("*.list"), path("${taxid}_sra.count")

    script:
    """
    bash search_sra.sh -t ${taxid} -n ${name} -k ${params.ncbi_api_key}
    """
}
