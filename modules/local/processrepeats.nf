process PROCESS_REPEATS {
    tag "PROCESS_REPEATS_${meta.id}"
    publishDir path: "${params.outdir}/results/${meta.id}/repetitive-elements", pattern: "*{masked,gff,tbl}"
    time '5d'
    label 'repeatmasker'

    input:
    tuple val(meta), path(masked)
    tuple val(meta), path(cat)
    tuple val(meta), path(library)

    output:
    tuple val(meta), path("${meta.id}.masked"),      emit: masked
    tuple val(meta), path("${meta.id}.gff"),         emit: gff
    tuple val(meta), path("${meta.id}.align"),       emit: align
    tuple val(meta), path("${meta.id}.cat"),         emit: cat
    tuple val(meta), path("${meta.id}.out"),         emit: out
    tuple val(meta), path("${meta.id}.tbl"),         emit: tbl

    script:
    """
    cat $cat > ${meta.id}.cat
    cp $masked ${meta.id}.masked

    ProcessRepeats -a -gff -lib $library -xsmall ${meta.id}.cat
    """
}
