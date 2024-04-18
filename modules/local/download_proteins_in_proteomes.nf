process DOWNLOAD_PROTEINS_IN_PROTEOMES {
    scratch true
    tag "DOWNLOAD_PROTEINS_IN_PROTEOMES_${name}"
    cache 'lenient'

    input:
    tuple val(name), val(meta), path(id_list)

    output:
    tuple val(name), val(meta), path{"${meta.taxid}_proteins.fa"}

    script:
    """
    # get all proteome IDs in the given clade
    curl --output ${meta.taxid}_UPIDs.list.zip \
        "https://rest.uniprot.org/proteomes/stream?compressed=true&format=list&query=%28%28taxonomy_id%3A${meta.taxid}%29%29+AND+%28proteome_type%3A1%29"

    mkdir -p proteomes

    for upid in \$(zcat ${meta.taxid}_UPIDs.list.zip); do
    # get all protein IDs in a proteome
    curl --output proteomes/\${upid}_protIDs.fa.zip \
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28proteome%3A\${upid}%29%29+AND+%28reviewed%3Atrue%29"
    done

    # merge all protein IDs
    zcat proteomes/* > ${meta.taxid}_proteins.fa
    """

    stub:
    """
    touch ${meta.taxid}_proteins.fa
    """
}
