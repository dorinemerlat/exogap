process GET_TAXONOMIC_LINEAGE {
    debug
    tag "GET_TAXONOMIC_LINEAGE_${meta.id}"


    conda (params.enable_conda ? 'bioconda requests==2.26.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/requests:2.26.0':
        'quay.io/biocontainers/requests:2.26.0' }"

    input:
    tuple val(meta), val(genome)

    output:
    stdout

    script:
    """
    get_taxonomic_lineage.py -t $meta.taxid -f $genome
    """
}
