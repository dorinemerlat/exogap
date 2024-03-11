process DOWNLOAD_PROTEINS {
    tag "DOWNLOAD_PROTEINS_${name}"
    publishDir "out/data/proteins/${taxid}/*", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(name), val(taxid), path(id_list)

    output:
    tuple val(taxid), val(name), path{"${taxid}_proteins.fa"}

    script:
    """
    curl --output ${taxid}_proteins.fa.zip \
        "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A${taxid}%29%29"

    if [[ \$(file ${taxid}_proteins.fa.zip) == *"gzip compressed"* ]]; then
        zcat ${taxid}_proteins.fa.zip > ${taxid}_proteins.fa
    else
        for i in \$(echo "1/100 101/200 201/300 301/400 401/600 601/800 801/*"); do
            interval=\$(echo \$i |sed "s@/@+TO+@g")
           curl --output ${taxid}_proteins_subset_\${interval}.fa.zip \
                "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A${taxid}%29%29+AND+%28length%3A%5B\${interval}%5D%29"
        done

        zcat ${taxid}_proteins_subset_*.fa.zip > ${taxid}_proteins.fa
    fi
    """

    stub:
    """
    touch ${taxid}_proteins.fa
    """
}

