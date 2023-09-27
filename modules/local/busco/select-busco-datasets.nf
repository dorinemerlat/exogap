process SELECT_BUSCO_DATASETS {
    tag "SELECT_BUSCO_DATASETS_${genome_id}"

    conda (params.enable_conda ? 'bioconda busco==5.5.0--pyhdfd78af_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.0--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.5.0--pyhdfd78af_0' }"

    input:
    path busco_dataset
    tuple val(genome_id), path(genome_lineage)

    output:
    stdout

    script:
    """
    awk 'NR==2' $genome_lineage | tr "\t" "\n"  |grep -if - $busco_dataset |tr "\n" "," |sed "s/,\$//g"
    """
}
