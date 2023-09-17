process REPEATMODELER {
    tag "REPEATMODELER_$genome_id"
    cpus 32
    time '20d'

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.3--pl5321h9ee0642_0':
        'quay.io/biocontainers/repeatmodeler:2.0.3--pl5321h9ee0642_0' }"

    input:
    tuple val(genome_id), path(genome_path)

    output:
    tuple val("${genome_id}-families"), path("${genome_id}-families.fa")

    script:
    """
    BuildDatabase -name $genome_id -engine ncbi $genome_path

    RepeatModeler -pa $task.cpus -engine ncbi -database $genome_id -LTRStruct
    """
}
