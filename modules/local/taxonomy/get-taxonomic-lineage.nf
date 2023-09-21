process GET_TAXONOMIC_LINEAGE {
    tag "GET_TAXONOMIC_LINEAGE_${genome_id}"

    // conda (params.enable_conda ? 'bioconda requests==2.26.0' : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/requests:2.26.0':
    //     'quay.io/biocontainers/requests:2.26.0' }"

    input:
    tuple val(genome_id), val(taxid)

    output:
    tuple val(genome_id), path("${genome_id}_lineage.tsv")

    script:
    """
    get_taxonomic_lineage.py -t $taxid -n $genome_id -o ${genome_id}_lineage.tsv
    """
}
