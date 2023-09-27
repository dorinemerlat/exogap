process BUSCO_SET_SELECTION {
    tag "BUSCO_SET_SELECTION_${genome_id}"
    debug true

    conda (params.enable_conda ? 'bioconda busco==5.5.0--pyhdfd78af_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.0--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.5.0--pyhdfd78af_0' }"

    input:
    path (dataset_list)
    tuple val(genome_id), path(genome_path), val(taxid)

    output:
    // tuple val(genome_id), path(genome_path),    emit: genome
    tuple val(genome_id), stdout

    script:
    """
    curl -X GET "https://lbgi.fr/api/taxonomy/lineage/${taxid}" -H  "accept: application/json" | python -m json.tool > ${genome_id}_lineage.txt

    grep -if list_datasets_reformated.txt ${genome_id}_lineage.txt |awk -F '"' '{print \$(NF-1)}' |grep -if - list_datasets.txt |awk '{print \$NF}'
    """
}
