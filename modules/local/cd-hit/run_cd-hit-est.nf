process RUN_CD_HIT_EST {
    tag "RUN_CD_HIT_EST_$input"
    cpus 10
    memory '16 GB'
    time '1d'

    conda (params.enable_conda ? 'cd-hit==4.8.1--h5b5514e_8' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5b5514e_7':
        'biocontainers/cd-hit:v4.6.8-2-deb_cv1' }"

    input:
    tuple val(input), val(output)

    output:
    tuple val(input), path("${output}")

    script:
    """
    memory=\$(echo $task.memory | sed 's/Gb//')
    memory=\$((\$memory * 1024))

    cd-hit-est -i $input -o $output -c 0.95 -g 1 -n 8 -d 10 \$memory -T $task.cpus
    """
}
