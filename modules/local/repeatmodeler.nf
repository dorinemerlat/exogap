process REPEATMODELER {
    scratch true
    tag "REPEATMODELER_${id}"
    cpus 30
    time '20d'
    label 'repeatmodeler'

    input:
    tuple  val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}_families.fa")

    script:
    """
    BuildDatabase -name $id -engine ncbi $genome

    RepeatModeler -pa $task.cpus -engine ncbi -database $id -LTRStruct
    """

    stub:
    """
    touch ${id}_families.fa
    """
}
