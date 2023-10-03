process PROCESS_REPEATS {
    tag "PROCESS_REPEATS_${meta.id}"
    time '1d'
    publishDir path: "results/${meta.id}/repeatmasker", pattern: "*{masked,gff}"

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'quay.io/biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(masked)
    tuple val(meta), path(cat)
    tuple val(library_id), path(library_path)

    output:
    tuple val(meta), path("${meta.id}.masked"),      emit: masked
    tuple val(meta), path("${meta.id}.gff"),         emit: gff
    tuple val(meta), path("${meta.id}-complex.gff"), emit: gff_complex
    tuple val(meta), path("${meta.id}.align"),       emit: align
    tuple val(meta), path("${meta.id}.cat"),         emit: cat

    script:
    """
    cat $cat > ${meta.id}.cat
    cp $masked ${meta.id}.masked

    ProcessRepeats -a -gff -lib $repeats_lib -xsmall ${meta.id}.cat

    # Isolation of complex repeats
    grep -v -e "Satellite" -e ")n" -e "-rich"  ${meta.id}.gff > ${meta.id}-complex.gff.unreformated

    # Reformatting to make it work with Maker
    cat ${meta.id}-complex.gff.unreformated | \
        perl -ane '\$id; if(!/^\#/){{@F = split(/\t/, \$_); chomp \$F[-1];\$id++; \$F[-1] .= "\;ID=\$id"; \$_ =join("\t", @F)."\n"}} print \$_' \
        > ${meta.id}-complex.gff
    """
}
