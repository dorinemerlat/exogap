process REPEATMASKER {
    tag "REPEATMASKER_${meta.id}"
    cpus 30
    time '20d'
    label 'repeatmasker'

    input:
    tuple val(meta), path(genome)
    tuple val(meta), path(library)

    output:
    tuple val(meta), path("${genome}.masked"),    emit: masked
    tuple val(meta), path("${genome}.cat"),       emit: cat

    script:
    """
    RepeatMasker -pa $task.cpus -e ncbi -nolow -lib $library -a -gccalc -norna -excln \
        -gff -s -xsmall -gccalc -excln -gff -s -html $genome

    gunzip *.gz
    """
}
