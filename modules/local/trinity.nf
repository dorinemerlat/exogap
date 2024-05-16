process TRINITY {
    tag "TRINITY_${specie_name}"
    cpus 40
    memory "20 GB"
    time '5d'
    label 'trinity'

    input:
    tuple val(clade_taxid), val(specie_name), val(specie_taxid), path(fastq1), path(fastq2)

    output:
    tuple val(clade_taxid), val(specie_name), val(specie_taxid), path("trinity_${specie_name}.Trinity.fasta")

    script:
    """
    replace_spaces_with_commas() {
        echo \$1 | sed 's/ /,/g'
    }

    fastq1_comma=\$(replace_spaces_with_commas $fastq1)
    fastq2_comma=\$(replace_spaces_with_commas $fastq2)
    memory=\$(echo $task.memory | sed 's/ GB/G/')

    Trinity --seqType fq --left \${fastq1_comma} --right \${fastq2_comma} --output trinity_${specie_name} --CPU $task.cpus --max_memory \$memory
    """

    stub:
    """
    touch trinity_${specie_name}.Trinity.fasta
    """
}
