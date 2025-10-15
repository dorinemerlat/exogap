process HITE {
    tag "${id}"
    cpus 10
    time '10d'
    label 'hite'
    scratch '/data/merlat'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("out")

    script:
    """
    python /HiTE/main.py \
        --genome ${genome} \
        --thread 60 \
        --out_dir out \
        --domain 1 \
        --plant 0 \
        --recover 1 \
        --search_struct 1
    """

    stub:
    """
    touch ${id}_earlgrey_families.fa
    """
}
