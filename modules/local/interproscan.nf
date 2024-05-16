process INTERPROSCAN {
    tag "INTERPROSCAN_${id}"
    label 'interproscan'
    cpus 20

    input:
    tuple val(id), val(meta), path(proteins)

    output:
    tuple val(id), val(meta), path("${id}.ips.xml"),    emit: xml
    tuple val(id), val(meta), path("${id}.ips.tsv"),    emit: tsv
    tuple val(id), val(meta), path("${id}.ips.gff3"),   emit: gff

    script:
    """
    interproscan.sh -i ${proteins} -b ${id}.ips -f TSV,XML,GFF3 -cpu 1 -dp -goterms -iprlookup -pa -t p --cpu ${task.cpus}
    """

    stub:
    """
    touch ${id}.ips.xml ${id}.ips.tsv ${id}.ips.gff3
    """
}
