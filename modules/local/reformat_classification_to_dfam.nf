process REFORMAT_CLASSIFICATION_TO_DFAM {
    scratch true
    tag "REFORMAT_CLASSIFICATION_TO_DFAM_${id}"

    input:
    tuple  val(id), val(meta), path(gff)
    tuple  val(id), val(meta), path(cat)
    tuple  val(id), val(meta), path(out)
    tuple  val(id), val(meta), path(tbl)
    tuple  val(id), val(meta), path(align)
    each path(dfam)

    output:
    tuple  val(id), val(meta), path("${id}_reformat.gff"),      emit: gff
    tuple  val(id), val(meta), path("${id}_reformat.cat"),      emit: cat
    tuple  val(id), val(meta), path("${id}_reformat.out"),      emit: out
    tuple  val(id), val(meta), path("${id}_reformat.tbl"),      emit: tbl
    tuple  val(id), val(meta), path("${id}_reformat.align"),    emit: align

    script:
    """
    touch ${id}_reformat.gff ${id}_reformat.cat ${id}_reformat.out ${id}_reformat.tbl ${id}_reformat.align
    """

    stub:
    """
    touch ${id}_reformat.gff ${id}_reformat.cat ${id}_reformat.out ${id}_reformat.tbl ${id}_reformat.align
    """
}
