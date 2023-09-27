process REPEATMASKER {
    tag "REPEATMASKER_$genome_id"
    cpus 30
    time '20d'
    publishDir "results/$genome_id/repeatmasker/$library_id"

    conda (params.enable_conda ? 'repeatmasker==4.1.5--pl5321hdfd78af_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'quay.io/biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(genome_id), path(genome_path)
    tuple val(library_id), path(library_path)

    output:
    tuple val(genome_id), path("${genome_path}.masked"),    emit: masked
    path "${genome_path}.cat",                              emit: cat

    script:
    """
    RepeatMasker -pa $task.cpus -e ncbi -nolow -lib $library_path -a -gccalc -norna -excln \
        -gff -s -xsmall -gccalc -excln -gff -s -html $genome_path

    gunzip *.gz
    """
}
