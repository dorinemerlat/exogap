process CALCULATE_GENOME_SIZE {
    tag "$id"
    cache 'lenient'
    label 'bioawk'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta) , path("${id}.size"),  emit: size
    path "${id}_size.csv",                          emit: size_for_info

    script:
    """
    bioawk -c fastx '{print length(\$seq)}' $genome |awk -v id=$id 'BEGIN{OFS=",";} {sum+=\$1} END {print id, sum;}' > "${id}.size"
    awk 'BEGIN {print "id,genome_size"} {print}' ${id}.size > ${id}_size.csv
    """

    stub:
    """
    echo "${id}, 0" > ${id}.size
    awk 'BEGIN {print "id,genome_size"} {print}' ${id}.size > ${id}_size.csv
    """
}
