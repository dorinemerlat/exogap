process GENOME_FORMAT {
    tag "GENOME_FORMAT_$genome_id"
    publishDir "results/$genome_id"

    input:
    tuple val(genome_id), path(genome_path)

    output:
    tuple val(genome_id), path(genome_path),            emit: fasta
    tuple val(genome_id), path("${genome_id}.json"),    emit: json

    script:
    """
    mv $genome_path ${genome_id}_unreformated.fa
    reformat_genome.py -i ${genome_id}_unreformated.fa -o ${genome_id} -p ${genome_id}-
    """
}
