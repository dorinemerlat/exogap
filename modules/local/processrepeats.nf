process PROCESS_REPEATS {
    scratch true
    tag "PROCESS_REPEATS_${id}"
    time '5d'
    label 'repeatmasker'

    input:
    tuple  val(id), val(meta), path(masked)
    tuple  val(id), val(meta), path(cat)
    each path(library)

    output:
    tuple  val(id), val(meta), path("${id}.masked"),      emit: masked
    tuple  val(id), val(meta), path("${id}.gff"),         emit: gff
    tuple  val(id), val(meta), path("${id}.align"),       emit: align
    tuple  val(id), val(meta), path("${id}.cat"),         emit: cat
    tuple  val(id), val(meta), path("${id}.out"),         emit: out
    tuple  val(id), val(meta), path("${id}.tbl"),         emit: tbl

    script:
    """
    cat $cat > ${id}.cat
    cp $masked ${id}.masked

    ProcessRepeats -a -gff -lib $library -xsmall ${id}.cat
    """

    stub:
    """
    touch ${id}.masked ${id}.gff ${id}.align ${id}.cat ${id}.out ${id}.tbl
    """
}
