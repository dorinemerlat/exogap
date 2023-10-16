process DOWNLOAD_BUSCO_DATASETS {
    tag "DOWNLOAD_BUSCO_DATASETS"
    publishDir "${params.outdir}/results/programs-outputs/busco/", mode: 'symlink'

    conda (params.enable_conda ? 'bioconda busco==5.5.0--pyhdfd78af_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.0--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.5.0--pyhdfd78af_0' }"

    output:
    path "BUSCO_datasets.txt"

    script:
    """
    busco --list-datasets | grep " - " |  awk '{print \$NF}' > BUSCO_datasets.txt
    """
}
