process GET_TAXONOMIC_LINEAGE {
    tag "GET_TAXONOMIC_LINEAGE_${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://dorinemerlat/python-exogap:v1.02':
    'dorinemerlat/python-exogap:v1.02' }"

    input:
    tuple val(meta), val(genome)

    output:
    stdout

    script:
    """
    get_taxonomic_lineage.py -t $meta.taxid -f $genome
    """
}
