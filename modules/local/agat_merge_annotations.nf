process AGAT_MERGE_ANNOTATIONS {

    tag "AGAT_MERGE_ANNOTATIONS_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff_files), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("merged.gff"), val(iteration)

    script:
    """
    gff_args=""
    for gff in $gff_files; do
    gff_arg=\$(echo "\$gff_args --gff \$gff")
    done

    agat_sp_keep_longest_isoform.pl \$gff_args -o merged.gff
    """

    stub:
    """
    touch merged.gff
    """
}
