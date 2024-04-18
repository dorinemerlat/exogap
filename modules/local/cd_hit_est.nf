process CD_HIT_EST {
    scratch true
    tag "CD_HIT_EST_$library_id"
    cpus 64
    time '1d'
    label 'cd_hit'

    input:
    tuple val(library_id), val(meta), path(library_path), val(word_length), val(lenght_of_throw_away_sequence), val(alignment_coverage_for_the_shorter_sequence)

    output:
    tuple val("${library_id}"), val(meta), path("${library_path}.nr")

    script:
    """
    memory=\$(echo "$task.memory" | sed 's/GB//')
    memory=\$((\$memory * 1024))

    cd-hit-est -i $library_path -o ${library_path}.nr \
        -c 0.95 \
        -g 1 \
        -n $word_length \
        -l $lenght_of_throw_away_sequence \
        -M \$memory \
        -T $task.cpus \
        -aS $alignment_coverage_for_the_shorter_sequence
    """

    stub:
    """
    touch ${library_path}.nr
    """
}

    // cd-hit-est -i $library_path -o ${library_path}.nr \
    //     -c 0.95 \ # sequence identity threshold
    //     -g 1 \ # if set to 1, the program will cluster it into the most similar cluster that meet the threshold
    //     -n $word_length \ #10 for repetitive elements, 5 for others
    //     -l $lenght_of_throw_away_sequence \ # 80 for repetitive elements, 10 for others
    //     -M \$memory \
    //     -T $task.cpus \
    //     -aS $alignment_coverage_for_the_shorter_sequence \ # 98 for repetitive elements, 0 for others
