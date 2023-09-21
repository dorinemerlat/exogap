process BUSCO {
    tag "BUSCO_${genome_id}_${dataset}_${mode}"
    cpus 15
    time '1d'
    publishDir "results/$genome_id/genome_quality"

    conda (params.enable_conda ? 'bioconda busco==5.5.0--pyhdfd78af_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.0--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.5.0--pyhdfd78af_0' }"

    input:
    tuple val(genome_id), path(genome_path)
    each dataset
    val mode

    output:
    path "${genome_id}_${dataset}.json"

    script:
    """
    echo busco -i $genome_path -l $dataset -o "${genome_id}_${dataset}.json" -m $mode -c N $task.cpus
    """
}
