process TRNASCANSE {
    tag "${id}"
    label 'trnascan'
    cpus 20

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}_rnammer.gff"), emit: gff
    tuple val(id), val(meta), path("${id}_rnammer.fa"), emit: fasta
    tuple val(id), val(meta), path("${id}_rnammer.out"), emit: out
    tuple val(id), val(meta), path("${id}_rnammer.struct"), emit: struct
    tuple val(id), val(meta), path("${id}_rnammer.stats"), emit: stats

    script:
    """
    tRNAscan-SE -E \
        --thread ${task.cpus} \
        --output ${id}_rnammer.out \
        --struct ${id}_rnammer.struct \
        --gff ${id}_rnammer.gff \
        --fasta ${id}_rnammer.fa \
        --stats ${id}_rnammer.stats \
        --progress \
        ${genome}
    """

    stub:
    """
    touch ${id}_rnammer.gff ${id}_rnammer.fa ${id}_rnammer.out ${id}_rnammer.struct ${id}_rnammer.stats
    """
}
