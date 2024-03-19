process MAKER_BY_SIMILARITY {
    tag "MAKER_BY_SIMILARITY_${genome}_${type}"
    label 'maker'
    time '20d'

    input:
    tuple val(id), val(meta), path(genome), path(ctl), val(type)

    output:
    tuple  val(id), val(meta), path("${id}/${id}.maker.output/${id}_master_datastore_index.log")

    script:
    """
    mpiexec -wdir . -n $task.cpus maker -base $id -fix_nucleotides
    """

    stub:
    """
    touch maker_opts.ctl
    """
}
