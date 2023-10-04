process REPEATMASKER {
    tag "REPEATMASKER_${meta.id}"
    cpus 30
    time '20d'

    conda (params.enable_conda ? 'repeatmasker==4.1.5--pl5321hdfd78af_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'quay.io/biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

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
