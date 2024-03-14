process DOWNLOAD_PROTEINS {
    tag "DOWNLOAD_PROTEINS_${name}"
    cache 'lenient'

    input:
    tuple val(name), val(meta), path(id_list)

    output:
    tuple val(name), val(meta), path{"${meta.taxid}_proteins.fa"}

    script:
    """
    curl --output ${meta.taxid}_proteins.fa.zip \
        "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A${meta.taxid}%29%29"

    if [[ \$(file ${meta.taxid}_proteins.fa.zip) == *"gzip compressed"* ]]; then
        zcat ${meta.taxid}_proteins.fa.zip > ${meta.taxid}_proteins.fa
    else
        for i in \$(echo "1/100 101/200 201/300 301/400 401/600 601/800 801/*"); do
            interval=\$(echo \$i |sed "s@/@+TO+@g")
           curl --output ${meta.taxid}_proteins_subset_\${interval}.fa.zip \
                "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A${meta.taxid}%29%29+AND+%28length%3A%5B\${interval}%5D%29"
        done

        zcat ${meta.taxid}_proteins_subset_*.fa.zip > ${meta.taxid}_proteins.fa
    fi
    """

    stub:
    """
    touch ${meta.taxid}_proteins.fa
    """
}

