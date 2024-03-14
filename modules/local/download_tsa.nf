process DOWNLOAD_TSA {
    tag "DOWNLOAD_TSA_${name}"
    cache 'lenient'
    label 'sratools'

    input:
    tuple val(name), val(meta), path(id_list)

    output:
    tuple val(name), val(meta), path("${meta.taxid}_transcriptomes.fasta")

    script:
    """
    prefetch --option-file $id_list

    srr_directories=\$(find . -type f ! -name '.*')
    for srr in \${srr_directories} ; do
        fasterq-dump --fasta \${srr}
    done

    cat *fasta > ${meta.taxid}_transcriptomes.fasta
    """

    stub:
    """
    touch ${meta.taxid}_transcriptomes.fasta
    """
}
