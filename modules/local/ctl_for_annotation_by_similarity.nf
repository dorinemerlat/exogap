process CTL_FOR_ANNOTATION_BY_SIMILARITY {
    tag "CTL_FOR_ANNOTATION_BY_SIMILARITY_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(proteins), path(transcripts), path(repeats), val(iteration)

    output:
    tuple  val(id), val(meta), path(genome), path("maker_opts.ctl"), val(iteration)

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
    """

    stub:
    """
    touch maker_opts.ctl
    """
}
