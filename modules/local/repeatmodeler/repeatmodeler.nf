process REPEATMODELER {
    tag "REPEATMODELER_${meta.id}"
    cpus 30
    time '20d'

    conda (params.enable_conda ? 'repeatmodeler==2.0.3--pl5321h9ee0642_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.3--pl5321h9ee0642_0':
        'quay.io/biocontainers/repeatmodeler:2.0.3--pl5321h9ee0642_0' }"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}-families.fa")

    script:
    """
    BuildDatabase -name $meta.id -engine ncbi $genome

    RepeatModeler -pa $task.cpus -engine ncbi -database $meta.id -LTRStruct
    """
}
