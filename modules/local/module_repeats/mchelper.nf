process MCHELPER {
    tag "${id}"
    cpus 50
    time '20d'
    label 'mchelper'
    errorStrategy 'ignore'

    input:
    tuple val(id), val(meta), path(families), path(genome)

    output:
    tuple val(id), val(meta), path("classifiedModule/kept_seqs_classified_module.fa"), emit: classified
    tuple val(id), val(meta), path("unclassifiedModule/kept_seqs_unclassified_module.fa"), emit: unclassified

    script:
    """
    merged_input="${id}_merged_input.fa"
    cat $families > \$merged_input
    input=\$merged_input

    python3 MCHelper/MCHelper.py \
        -g ${genome} \
        -l \${input} \
        -o out_${id} \
        --input_type fasta \
        -r A \
        -t ${task.cpus} \
        -b arthropoda_odb10.hmm \
        -a F
    """

    stub:
    """
    touch classifiedModule/kept_seqs_classified_module.fa unclassifiedModule/kept_seqs_unclassified_module.fa
    """
}
