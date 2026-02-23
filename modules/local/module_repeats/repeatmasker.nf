process REPEATMASKER {
    tag "$id"
    cpus 30
    time '20d'
    // label 'repeatmasker'
    errorStrategy 'ignore'

    input:
    tuple val(id), val(meta), path(genome), path(library)

    output:
    tuple  val(id), val(meta), path("${genome}.masked"),    emit: masked
    // tuple  val(id), val(meta), path("${genome}.cat"),       emit: cat
    tuple  val(id), val(meta), path("${genome}.out"),       emit: out
    tuple  val(id), val(meta), path("${genome}.out.gff"),   emit: gff
    tuple  val(id), val(meta), path("${genome}.tbl"),       emit: tbl
    tuple  val(id), val(meta), path("${genome}.align"),     emit: align

    script:
    """
    RepeatMasker -pa $task.cpus -e ncbi -lib $library -a -gccalc -norna -excln  -s -xsmall -gff $genome

    gunzip *.gz
    """

    stub:
    """
    touch ${genome}.masked ${genome}.cat ${genome}.out ${genome}.out.gff ${genome}.tbl ${genome}.align
    """
}
