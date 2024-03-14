process CD_HIT {
    tag "CD_HIT_$library_id"
    cpus 10
    memory '16 GB'
    time '1d'
    label 'cd_hit'

    input:
    tuple val(library_id), val(meta), path(library_path)

    output:
    tuple val(library_id), val(meta), path("${library_path}.nr")

    script:
    """
    memory=\$(echo $task.memory | sed 's/Gb//')
    memory=\$((\$memory * 1024))

    cd-hit -i $library_id -o ${library_path}.nr -c 0.95 -n 5 -g 1 -d 0 -M \$memory -T $task.cpus
    """

    stub:
    """
    touch ${library_path}.nr
    """
}
