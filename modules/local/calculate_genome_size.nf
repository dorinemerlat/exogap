process calculate_genome_size {
    tag "calculate_genome_size_$id"
    cache 'lenient'
    label 'bioawk'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path(genome), path("${id}.csv")

    script:
    """
    bioawk -c fastx '{print length(\$seq)}' $genome |awk -v id=$id 'BEGIN{OFS=",";} {sum+=\$1} END {print id, sum;}' > "${id}.csv"
    """

    stub:
    """
    echo "${id}, 0;" > ${id}.csv
    """
}
