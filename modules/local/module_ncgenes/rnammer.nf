process RNAMMER {
    tag "RNAMMER_${id}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}_rnammer.gff")

    script:
    """
    /tempor/merlat/exogaptwo/bin/run_rnammer.sh ${id} ${genome}
    """

    stub:
    """
    touch ${id}_rnammer.gff ${id}_rnammer.fa
    """
}
