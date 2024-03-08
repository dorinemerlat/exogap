process CD_HIT {
    tag "CD_HIT_$input"
    cpus 10
    memory '16 GB'
    time '1d'
    label 'cd_hit'

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
