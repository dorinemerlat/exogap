process CD_HIT_FOR_REPEATS {
    tag "CD_HIT_FOR_REPEATS_$library_id"
    cpus 10
    time '1d'

    conda (params.enable_conda ? 'cd-hit==4.8.1--h5b5514e_8' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5b5514e_7':
        'biocontainers/cd-hit:v4.6.8-2-deb_cv1' }"

    input:
    tuple val(library_id), path(library_path)

    output:
    tuple val("${library_id}"), path("${library_id}-noredundance.fa")

    script:
    """
    cd-hit-est -i $library_path -o ${library_id}-noredundance.fa -c 0.95 -n 10 -g 1 -l 80 -aS 98 -d 10 -M 16000 -T 8
    """
}
