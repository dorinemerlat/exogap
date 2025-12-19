process BARRNAP {
    tag "${id}"
    // label 'barrnap'
    cpus 30

    input:
    tuple val(id), val(meta), path(genome), val(type)

    output:
    tuple val(id), val(meta), path("${id}_barrnap_${type}.gff"), emit: gff
    tuple val(id), val(meta), path("${id}_barrnap_${type}.fa"), emit: fasta
    script:
    """
    module load genomics/bedtools/
    /tempor/merlat/exogaptwo/bin/run_barrnap.sh ${id} ${genome} ${type} ${task.cpus}
    """

    stub:
    """
    touch ${id}_barrnap_${type}.gff
    """
}
