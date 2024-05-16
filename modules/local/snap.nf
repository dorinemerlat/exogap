process SNAP {
    tag "SNAP_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple  val(id), val(meta), path(genome), path("${id}_snap.hmm"), val(iteration)

    script:
    """
    maker2zff -x 0.25 -l 50 $gff

    # Gather some statistics and validate
    fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
    fathom genome.ann genome.dna -validate > validate.log 2>&1

    # Collects training sequences and annotations (over 1000 around) for training
    fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
    fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

    # Create training parameters
    mkdir params
    cd params
    forge ../export.ann ../export.dna > ../forge.log 2>&1
    cd ..

    # Assembling the HMM
    hmm-assembler.pl genome params > ${id}_snap.hmm
    """

    stub:
    """
    touch ${id}_snap.hmm
    """
}
