process plot_busco_summary {
    publishDir "${params.outdir}/results/programs-outputs/busco/", mode: 'symlink'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://r-base:4.3.1':
    'r-base:4.3.1' }"

    input:
    tuple val(ids), val(metas), path(newick), path(csv)

    output:
    tuple val(name), path(data),  emit: data
    tuple val(name), path("${name}.png"),  emit: plots

    script:
    """
    R $data
    """
}
