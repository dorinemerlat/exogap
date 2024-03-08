process BUSCO {
    tag "BUSCO_${id}_${dataset}_${mode}"
    cpus 15
    time '1d'
    publishDir "${params.outdir}/results/programs-outputs/busco/", mode: 'symlink'
    label: 'busco'

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

    stub:
    """
    touch short_summary.*.json
    touch short_summary.*.txt
    """
}
