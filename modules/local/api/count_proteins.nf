process COUNT_PROTEINS {
    tag "COUNT_PROTEINS_${name}"
    publishDir "out/data/proteins/${taxid}/*", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(taxid), val(name)

    output:
    tuple val(taxid), val(name), path{"${taxid}_proteins.list.zip"}, path{"${taxid}_proteins_count.list"}, path{"${taxid}_proteins_empty.txt"}

    script:
    """
    curl --output ${taxid}_proteins.list.zip \
        "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=list&query=%28%28taxonomy_id%3A${taxid}%29%29"

    if [[ \$(file ${taxid}_proteins.list.zip) == *"gzip compressed"* ]]; then
        zcat ${taxid}_proteins.list.zip > ${taxid}_proteins.list
    else
        for i in \$(echo "1/100 101/200 201/300 301/400 401/600 601/800 801/*"); do
            interval=\$(echo \$i |sed "s@/@+TO+@g")
           curl --output ${taxid}_proteins_subset_\${interval}.list.zip \
                "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=list&query=%28%28taxonomy_id%3A${taxid}%29%29+AND+%28length%3A%5B\${interval}%5D%29"
        done

        zcat ${taxid}_proteins_subset_*.list.zip > ${taxid}_proteins.list

        # if always too many results, don't use this parent taxid
        line_count=\$(cat ${taxid}_proteins.list.zip | grep -q "Too many results to retrieve" |wc -l)
        if [  \$line_count = 0 ]; then
            echo '' > ${taxid}_proteins.list
        fi
    fi

    # zcat ${taxid}_proteins.list.zip > ${taxid}_proteins.list

    proteins_count=\$(wc -l < ${taxid}_proteins.list)
    echo "${taxid},${name},\${proteins_count}" > ${taxid}_proteins_count.list
    touch "${taxid}_proteins_empty.txt"
    """

    // # stub:
    // # """
    // # touch ${taxid}_proteins.list.zip
    // # echo "${taxid},${name},0" > ${taxid}_proteins_count.list
    // # """
}

