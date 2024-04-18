process CREATE_AUGUSTUS_SET {
    scratch true
    tag "CREATE_AUGUSTUS_SET_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("genes.gb.train"), path("genes.gb.test"), val(iteration)

    script:
    """
    gff2gbSmallDNA.pl $gff $genome 1000 genes.gb |tee gff2gbSmallDNA.log
    randomSplit.pl genes.gb 100
    """

    stub:
    """
    touch genes.gb.train genes.gb.test
    """
}
