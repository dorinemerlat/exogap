process REPEATMASKER {
    tag "REPEATMASKER_$genome_id"
    cpus 30
    time '20d'

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'quay.io/biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(genome_id), path(genome_path)
    path library_path

    output:
    tuple val(genome_id), path("${genome_path}.masked"),    emit: masked
    tuple val(genome_id), path("${genome_path}.tbl"),       emit: tbl
    tuple val(genome_id), path("${genome_path}.out"),       emit: out
    tuple val(genome_id), path("${genome_path}.align"),     emit: align
    tuple val(genome_id), path("${genome_path}.cat"),       emit: cat

    script:
    """
    RepeatMasker -pa $task.cpus -e ncbi -nolow -lib $library_path -a -gccalc -norna -excln \
        -gff -s -xsmall -gccalc -excln -gff -s -html $genome_path

    gunzip *.gz
    """
}
