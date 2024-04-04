process CREATE_AUGUSTUS_SET {
    tag "CREATE_AUGUSTUS_SET_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(gff)
    each val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("genes.gb.train"), path("genes.gb.test")

    script:
    """
    gff2gbSmallDNA.pl gff $genome 1000 genes.gb |tee gff2gbSmallDNA.log
    randomSplit.pl genes.gb 100

    """

    stub:
    """
    touch augustus.gff
    """
}-
