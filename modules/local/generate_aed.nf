process GENERATE_AED {
    scratch true
    tag "GENERATE_AED_${id}_${iteration}"
    label 'perl'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple  val(id), val(meta), path(genome), path("${id}_aed.stats"), val(iteration)

    script:
    """
    perl AED_cdf_generator.pl -b 0.020 $gff |grep '^[0,1]' > ${id}_aed.stats
    """

    stub:
    """
    touch ${id}_aed.stats
    """
}
