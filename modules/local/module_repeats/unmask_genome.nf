process UNMASK_GENOME {
    tag "${id}"

    input:
    tuple  val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path(genome)

    script:
    """
    mv $genome ${id}_softmasked.fa 
    awk '{if (\$0 ~ /^>/) print \$0; else print toupper(\$0)}' ${id}_softmasked.fa > $genome
    """

    stub:
    """
    """
}
