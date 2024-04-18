process MAKER_BY_SIMILARITY {
    scratch true
    tag "MAKER_BY_SIMILARITY_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(proteins), path(transcripts), path(repeats), val(iteration)

    output:
    tuple  val(id), val(meta), path(genome), path("${id}/${id}.maker.output/${id}_master_datastore_index.log"), val(iteration)

    script:
    """
    maker -CTL

    sed -i -e "s|^genome=|genome=${genome}|g" maker_opts.ctl
    sed -i -e "s|^est=|est=${transcripts}|g" maker_opts.ctl
    sed -i -e "s|^protein=|protein=${proteins}|g" maker_opts.ctl
    sed -i -e "s|^model_org=all|model_org=simple|g" maker_opts.ctl
    sed -i -e "s|^rm_gff=|rm_gff=${repeats}|g" maker_opts.ctl
    sed -i -e "s|^est2genome=0|est2genome=1|g" maker_opts.ctl
    sed -i -e "s|^protein2genome=0|protein2genome=1|g" maker_opts.ctl

    # run maker
    mpiexec -wdir . -n $task.cpus maker -base $id -fix_nucleotides
    """

    stub:
    """
    mkdir -p ${id}/${id}.maker.output
    touch ${id}/${id}.maker.output/${id}_master_datastore_index.log
    """
}
