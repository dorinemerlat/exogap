process BARRNAP {
    tag "BARRNAP_${id}_${type}"
    label 'barrnap'
    cpus 40

    input:
    tuple val(id), val(meta), path(genome), val(type)

    output:
    tuple val(id), val(meta), path("${id}_barrnap_${type}.gff")

    script:
    """
    barrnap --quiet --kingdom ${type} --threads $task.cpus --outseq ${id}_barrnap_${type}.fa $genome > ${id}_barrnap_${type}.gff.tmp

    grep -v '^#' ${id}_barrnap_${type}.gff.tmp \
        | sed -E 's/barrnap:[0-9]+(\.[0-9]+)?/barrnap/g' \
        | awk '{ nb++; print \$1"\t"\$2"\tgene\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\tID=barrnap_"nb";"\$9";family=Gene,rRNA\n" \
            \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\tID=rnammer_"nb"-rRNA;Parent=barrnap_"nb";"\$9";family=Gene,rRNA"}' \
        > ${id}_barrnap_${type}.gff
    """

    stub:
    """
    touch ${id}_barrnap_${type}.gff
    """
}
