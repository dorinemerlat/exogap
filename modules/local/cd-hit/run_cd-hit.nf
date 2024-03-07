process RUN_CD_HIT {
    tag "RUN_CD_HIT_$input"
    cpus 10
    memory '16 GB'
    time '1d'

    conda (params.enable_conda ? 'cd-hit==4.8.1--h5b5514e_8' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5b5514e_7':
        'biocontainers/cd-hit:v4.6.8-2-deb_cv1' }"

    input:
    tuple val(id), path(input), val(output)

    output:
    tuple val(id), path("${output}")

    script:
    """
    memory=\$(echo $task.memory | sed 's/Gb//')
    memory=\$((\$memory * 1024))

    cd-hit -i $input -o $output -c 0.95 -n 5 -g 1 -d 0 -M \$memory -T $task.cpus
    """
}
