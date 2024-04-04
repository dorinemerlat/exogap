process AGAT_FILTER_FEATURES_BY_ATTRIBUTE_VALUE {
    tag "AGAT_FILTER_FEATURES_BY_ATTRIBUTE_VALUE_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration), val(aed_threshold)

    output:
    tuple val(id), val(meta), path(genome), path("aed.gff"), val(iteration)

    script:
    """
    agat_sp_filter_feature_by_attribute_value.pl --gff $gff -o aed.gff --value $aed_threshold --a _AED -t ">="
    """

    stub:
    """
    touch aed.gff
    """
}
