process cd_hit_est {
    tag "cd_hit_est_$library_id"
    cpus 10
    time '1d'
    label 'cd_hit'

    input:
    tuple val(library_id), path(library_path), val(word_length), val(lenght_of_throw_away_sequence), val(alignment_coverage_for_the_shorter_sequence)

    output:
    tuple val("${library_id}"), path("${library_path}.nr")

    script:
    """
    cd-hit-est \
        -i $library_path \
        -o ${library_path}.nr \
        -c 0.95 \ # sequence identity threshold
        -g 1 \ # if set to 1, the program will cluster it into the most similar cluster that meet the threshold
        -n $word_length \ #10 for repetitive elements, 5 for others
        -l $lenght_of_throw_away_sequence \ # 80 for repetitive elements, 10 for others
        -M \$memory \
        -T $task.cpus \
        -aS $alignment_coverage_for_the_shorter_sequence \ # 98 for repetitive elements, 0 for others
    """
}
