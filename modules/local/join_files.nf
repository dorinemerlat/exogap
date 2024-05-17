process JOIN_FILES {
    // in tag, replace characters: [, ] by nothing
    tag "${id}"

    input:
    tuple val(id), path(file, stageAs: "file1.csv"), path(new_file), val(output_name)

    output:
    tuple val(id), path(output_name)

    script:
    """
    sort_csv_by_first_column() {
        { head -n 1 \$1 && tail -n +2 \$1 | sort -t ',' -k1,1; } > \$2
    }

    sort_csv_by_first_column file1.csv sorted_file1.csv
    sort_csv_by_first_column $new_file sorted_file2.csv

    join -t',' -1 1 -2 1 sorted_file1.csv sorted_file2.csv > $output_name
    """

    stub:
    """
    touch ${output_name}
    """
}