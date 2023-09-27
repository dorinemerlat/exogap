process RENAME_REPEATMODELER_OUTPUT {
    tag "RENAME_REPEATMODELER_OUTPUT_$genome_id"

        conda (params.enable_conda ? 'bioconda biopython==1.75' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(genome_id), path(genome_path)

    output:
    tuple val("${genome_id}"), path("${genome_id}-renamed.fa")

    script:
    """
    prefix=\$(echo $genome_id |sed -e "s/-families//g")
    rename_repeats.py -i $genome_path -o ${genome_id}-renamed.fa -p \$prefix
    """
}

