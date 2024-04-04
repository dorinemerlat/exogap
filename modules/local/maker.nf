process MAKER {
    tag "MAKER_${genome}_${iteration}"
    label 'maker'
    time '20d'

    input:
    tuple val(id), val(meta), path(genome), path(ctl), val(iteration)

    output:
    tuple  val(id), val(meta), path(genome), path("${id}/${id}.maker.output/${id}_master_datastore_index.log"), val(iteration)

    script:
    """
    mpiexec -wdir . -n $task.cpus maker -base $id -fix_nucleotides
    """

    stub:
    """
    mkdir -p ${id}/${id}.maker.output
    touch ${id}/${id}.maker.output/${id}_master_datastore_index.log
    """
}
