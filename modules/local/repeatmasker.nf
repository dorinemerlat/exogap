process REPEATMASKER {
    tag "REPEATMASKER_${id}"
    cpus 30
    time '20d'
    label 'repeatmasker'

    input:
    tuple  val(id), val(meta), path(genome)
    each path(library)

    output:
    tuple  val(id), val(meta), path("${genome}.masked"),    emit: masked
    tuple  val(id), val(meta), path("${genome}.cat"),       emit: cat

    script:
    """
    RepeatMasker -pa $task.cpus -e ncbi -nolow -lib $library -a -gccalc -norna -excln \
        -gff -s -xsmall -gccalc -excln -gff -s -html $genome

    gunzip *.gz
    """

    stub:
    """
    touch ${genome}.masked ${genome}.cat
    """
}
