process COUNT_PROTEINS_IN_PROTEOMES {
    tag "COUNT_PROTEINS_IN_PROTEOMES_${name}"

    publishDir "out/data/proteins_in_proteomes/${taxid}/*", mode: 'copy'

    cache 'lenient'

    input:
    tuple val(taxid), val(name)

    output:
    tuple val(taxid), val(name), path{"${taxid}_proteomes.list"}, path{"${taxid}_proteomes_count.list"}, path{"${taxid}_UPIDs.list.zip"}

    script:
    """
    # get all proteome IDs in the given clade
    curl --output ${taxid}_UPIDs.list.zip \
        "https://rest.uniprot.org/proteomes/stream?compressed=true&format=list&query=%28%28taxonomy_id%3A${taxid}%29%29+AND+%28proteome_type%3A1%29"

    line_count=\$(zcat ${taxid}_UPIDs.list.zip | wc -l)
    if [  \$line_count != 0 ]; then

        mkdir -p proteomes

        for upid in \$(zcat ${taxid}_UPIDs.list.zip); do
        # get all protein IDs in a proteome
            curl --output proteomes/\${upid}_protIDs.list.zip \
            "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=list&query=%28%28proteome%3A\${upid}%29%29+AND+%28reviewed%3Atrue%29" &
        done
        wait

        # merge all protein IDs
        zcat proteomes/* > ${taxid}_proteomes.list

        # count number of proteins for the given clade
        proteins_count=\$(wc -l < ${taxid}_proteomes.list)
        echo "${taxid},${name},\${proteins_count}" > ${taxid}_proteomes_count.list

    else

        touch ${taxid}_proteomes.list
        echo "${taxid},${name},0" > ${taxid}_proteomes_count.list
    fi
    """

    // stub:
    // """
    // echo "0" > ${taxid}_UPIDs.list && zip ${taxid}_UPIDs.list.zip ${taxid}_UPIDs.list
    // touch ${taxid}_proteomes.list
    // echo "${taxid},${name},0" > ${taxid}_proteomes_count.list
    // """
}
