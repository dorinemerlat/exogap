process CTL_FOR_AB_INITIO_ANNOTATION {
    tag "CTL_FOR_AB_INITIO_ANNOTATION_${genome}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(maker), path(snap), path(augustus)
    each val(iteration)s

    output:
    tuple  val(id), val(meta), path(genome), path("maker_opts.ctl"),      emit: masked


    script:
    """
    maker -CTL

    sed -i -e "s|^maker_gff=|maker_gff=${maker}|g" maker_opts.ctl
    sed -i -e "s|^est_pass=0|est_pass=1|g" maker_opts.ctl
    sed -i -e "s|^altest_pass=0|altest_pass=1|g" maker_opts.ctl
    sed -i -e "s|^protein_pass=0|protein_pass=1|g" maker_opts.ctl
    sed -i -e "s|^rm_pass=0|rm_pass=1|g" maker_opts.ctl
    sed -i -e "s|^model_pass=0|model_pass=1|g" maker_opts.ctl
    sed -i -e "s|^pred_pass=0|pred_pass=1|g" maker_opts.ctl
    sed -i -e "s|^other_pass=0|other_pass=1|g" maker_opts.ctl

    sed -i -e "s|^model_org=all|model_org=simple|g" maker_opts.ctl

    sed -i -e "s|^snaphmm=|snaphmm=${snap}|g" maker_opts.ctl
    sed -i -e "s|^augustus_species=|augustus_species=${augustus}|g" maker_opts.ctl

    sed -i -e "s|gmhmme3=|gmhmme3=${GeneMark}|g" maker_exe.ctl

    if [[ $Iteration == '3' ]] ; then # Itération = 1
        sed -i -e "s|^trna=0|trna=1|g" maker_opts.ctl
    fi
    """

    stub:
    """
    touch maker_opts.ctl
    """
}
