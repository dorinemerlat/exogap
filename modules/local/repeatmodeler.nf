process REPEATMODELER {
    tag "REPEATMODELER_${meta.id}"
    cpus 30
    time '20d'
    label 'repeatmodeler'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}-families.fa")

    script:
    """
    BuildDatabase -name $meta.id -engine ncbi $genome

    RepeatModeler -pa $task.cpus -engine ncbi -database $meta.id -LTRStruct
    """

    stub:
    """
    touch ${meta.id}-families.fa
    """
}
