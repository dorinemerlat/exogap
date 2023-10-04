process PROCESS_REPEATS {
    tag "PROCESS_REPEATS_${meta.id}"
    publishDir path: "${params.outdir}/results/${meta.id}/repetitive-elements", pattern: "*{masked,gff,tbl}"
    time '5d'

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'quay.io/biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(masked)
    tuple val(meta), path(cat)
    tuple val(meta), path(library)

    output:
    tuple val(meta), path("${meta.id}.masked"),      emit: masked
    tuple val(meta), path("${meta.id}.gff"),         emit: gff
    tuple val(meta), path("${meta.id}-complex.gff"), emit: gff_complex
    tuple val(meta), path("${meta.id}.align"),       emit: align
    tuple val(meta), path("${meta.id}.cat"),         emit: cat
    tuple val(meta), path("${meta.id}.out"),         emit: out
    tuple val(meta), path("${meta.id}.tbl"),         emit: tbl

    script:
    """
    cat $cat > ${meta.id}.cat
    cp $masked ${meta.id}.masked

    ProcessRepeats -a -gff -lib $library -xsmall ${meta.id}.cat

    # Isolation of complex repeats
    mv ${meta.id}.out.gff ${meta.id}.gff
    grep -v -e "Satellite" -e ")n" -e "-rich"  ${meta.id}.gff > ${meta.id}-complex.gff.unreformated

    # Reformatting to make it work with Maker
    reformat-repeats-gff.sh ${meta.id}-complex.gff.unreformated ${meta.id}-complex.gff
    """
}
