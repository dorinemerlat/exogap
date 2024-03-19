process SNAP {
    tag "SNAP${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(maker_gff), val(iteration)

    output:
    tuple  val(id), val(meta), path("${id}_snap.hmm"), val(iteration)

    script:
    """
    maker2zff -x 0.25 -l 50 $maker_gff

    # Rassembler quelques statistiques et valider
    fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
    fathom genome.ann genome.dna -validate > validate.log 2>&1

    # Collecte les séquences des training et les annotations (plus 1000 autour) pour le training
    fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
    fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

    # Créer les paramètres d'entrainement
    mkdir params
    cd params
    forge ../export.ann ../export.dna > ../forge.log 2>&1
    cd ..

    # Assembler le HMM
    hmm-assembler.pl genome params > ${id}_snap.hmm
    """

    stub:
    """
    touch ${id}_snap.hmm
    """
}
