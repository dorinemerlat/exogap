process MAKER_AB_INITIO {
    scratch true
    tag "MAKER_AB_INITIO_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(gff), path(snap), path(augustus), val(iteration)

    output:
    tuple  val(id), val(meta), path(genome), path("${id}/${id}.maker.output/${id}_master_datastore_index.log"), val(iteration)

    script:
    """
    # augustus config
    cp -r /usr/local/config .
    export AUGUSTUS_CONFIG_PATH=\$(pwd)/config
    cp $specie config/species/

    # maker config
    maker -CTL

    sed -i -e "s|^maker_gff=|maker_gff=${gff}|g" maker_opts.ctl
    sed -i -e "s|^est_pass=0|est_pass=1|g" maker_opts.ctl
    sed -i -e "s|^altest_pass=0|altest_pass=1|g" maker_opts.ctl
    sed -i -e "s|^protein_pass=0|protein_pass=1|g" maker_opts.ctl
    sed -i -e "s|^rm_pass=0|rm_pass=1|g" maker_opts.ctl
    sed -i -e "s|^model_pass=0|model_pass=1|g" maker_opts.ctl
    sed -i -e "s|^pred_pass=0|pred_pass=1|g" maker_opts.ctl
    sed -i -e "s|^other_pass=0|other_pass=1|g" maker_opts.ctl

    sed -i -e "s|^model_org=all|model_org=simple|g" maker_opts.ctl

    sed -i -e "s|^snaphmm=|snaphmm=${snap}|g" maker_opts.ctl
    sed -i -e "s|^augustus_species=|augustus_species=${id}|g" maker_opts.ctl

    sed -i -e "s|gmhmme3=|gmhmme3=${params.genemark}|g" maker_exe.ctl

    if [[ $iteration == '3' ]] ; then # Itération = 1
        sed -i -e "s|^trna=0|trna=1|g" maker_opts.ctl
    fi

    # run maker
    mpiexec -wdir . -n $task.cpus maker -base $id -fix_nucleotides
    """

    stub:
    """
    mkdir -p ${id}/${id}.maker.output
    touch ${id}/${id}.maker.output/${id}_master_datastore_index.log
    """
}
