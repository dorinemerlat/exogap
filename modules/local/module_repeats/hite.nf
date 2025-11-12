process HITE {
    tag "${id}"
    cpus 10
    time '10d'
    label 'hite'
    scratch '/data/merlat'
    errorStrategy 'ignore'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("out/confident_TE.cons.fa"), emit: confident_families
    tuple val(id), val(meta), path("out/low_confident_TE.cons.fa"), emit: low_confident_families
    tuple val(id), val(meta), path("out/confident_TE.cons.fa.domain"), emit: domains

    script:
    """
    python /HiTE/main.py \
        --genome ${genome} \
        --thread ${task.cpus} \
        --out_dir out \
        --domain 1 \
        --plant 0 \
        --recover 1 \
        --search_struct 1
    """

    stub:
    """
    touch out/confident_TE.cons.fa out/low_confident_TE.cons.fa out/confident_TE.cons.fa.domain
    """
}
