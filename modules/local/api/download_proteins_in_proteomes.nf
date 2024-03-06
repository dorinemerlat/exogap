process DOWNLOAD_PROTEINS_IN_PROTEOMES {
    tag "DOWNLOAD_PROTEINS_IN_PROTEOMES_${name}"

    publishDir "out/data/proteins_in_proteomes/${taxid}/*", mode: 'copy'

    cache 'lenient'

    input:
    tuple val(name), val(taxid), path(id_list)

    output:
    tuple val(taxid), val(name), path{"${taxid}_proteins.fa"}

    script:
    """
    # get all proteome IDs in the given clade
    curl --output ${taxid}_UPIDs.list.zip \
        "https://rest.uniprot.org/proteomes/stream?compressed=true&format=list&query=%28%28taxonomy_id%3A${taxid}%29%29+AND+%28proteome_type%3A1%29"

    mkdir -p proteomes

    for upid in \$(zcat ${taxid}_UPIDs.list.zip); do
    # get all protein IDs in a proteome
    curl --output proteomes/\${upid}_protIDs.fa.zip \
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28proteome%3A\${upid}%29%29+AND+%28reviewed%3Atrue%29"
    done

    # merge all protein IDs
    zcat proteomes/* > ${taxid}_proteins.fa
    """

    // stub:
    // """
    // echo "0" > ${taxid}_UPIDs.fa && zip ${taxid}_UPIDs.fa.zip ${taxid}_UPIDs.fa
    // touch ${taxid}_proteomes.fa
    // echo "${taxid},${name},0" > ${taxid}_proteomes_count.fa
    // """
}
