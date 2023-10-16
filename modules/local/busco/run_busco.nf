process BUSCO {
    tag "BUSCO_${id}_${dataset}_${mode}"
    cpus 15
    time '1d'
    publishDir "${params.outdir}/results/programs-outputs/busco/", mode: 'symlink'

    conda (params.enable_conda ? 'bioconda busco==5.5.0--pyhdfd78af_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.0--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.5.0--pyhdfd78af_0'}"

    input:
    tuple val(id), val(meta), path(genome), val(dataset), val(mode)

    output:
    tuple val(id), val(meta), path(genome), path("short_summary.*.json"), val(dataset), val(mode),  emit: json
    tuple val(id), val(meta), path(genome), path("short_summary.*.txt"),  val(dataset), val(mode),  emit: txt

    script:
    """
    busco -i $genome -l $dataset -o $id -m $mode -c $task.cpus
    mv ${id}/short_summary.* .
    """
}
