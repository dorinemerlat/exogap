process AGAT_MERGE_ANNOTATIONS {
    tag "AGAT_MERGE_ANNOTATIONS_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff1), path(gff2), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("merged.gff"), val(iteration)

    script:
    """
    agat_sp_keep_longest_isoform.pl --gff $gff1 --gff $gff2 -o merged.gff
    """

    stub:
    """
    touch merged.gff
    """
}
