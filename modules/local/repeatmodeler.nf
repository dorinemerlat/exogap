process REPEATMODELER {
    tag "REPEATMODELER_${id}"
    cpus 30
    time '20d'
    label 'repeatmodeler'

    input:
    tuple  val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}-families.fa")

    script:
    """
    BuildDatabase -name $id -engine ncbi $genome

    RepeatModeler -pa $task.cpus -engine ncbi -database $id -LTRStruct
    """

    stub:
    """
    touch ${id}-families.fa
    """
}
