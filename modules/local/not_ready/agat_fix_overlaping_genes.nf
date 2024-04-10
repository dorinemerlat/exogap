process AGAT_FIX_OVERLAPING_GENES {
    tag "AGAT_FIX_OVERLAPING_GENES_${id}_${iteration}"
    label 'agat'

    input:
    tuple val(id), val(meta), path(genome), path(gff), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path("distance.gff"), val(iteration)

    script:
    """
    agat_sp_fix_overlaping_genes.pl --gff $gff -o distance.gff -d 1000
    """

    stub:
    """
    touch distance.gff
    """
}
