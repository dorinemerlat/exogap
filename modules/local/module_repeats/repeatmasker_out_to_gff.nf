process REPEATMASKER_OUT_TO_GFF {
    tag "${id}"

    input:
    tuple val(id), val(meta), path(out), path(table)

    output:
    tuple val(id), path("${id}_repeats.gff")

    script:
    """
    repeatmasker_out_to_gff.py -r ${out} -m ${table} -o ${id}_repeats.gff
    """

    stub:
    """
    touch ${id}_repeats.gff
    """
}
