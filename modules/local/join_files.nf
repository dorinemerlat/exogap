process JOIN_FILES {
    tag "${id}"

    input:
    tuple val(id), path(files_to_join), val(output_name)

    output:
    tuple val(id), path(output_name)

    script:
    """
    # Define a function to sort CSV files by their first column
    sort_csv_by_first_column() {
        { head -n 1 \$1 && tail -n +2 \$1 | sort -t ',' -k1,1; } > \$2
    }

    # Function to rename 'info.csv' to 'info.csv.tmp' and return the file name
    rename_if_same_names() {
        local file_to_check=\$1

        if [ "\$file_to_check" == $output_name ]; then
            mv "\$file_to_check" "\${file_to_check}.tmp"
            echo "\${file_to_check}.tmp"
        else
            echo "\$file_to_check"
        fi
    }

    # Sort the first file and create the output file
    first_file=\$(echo $files_to_join | cut -d' ' -f1)
    sort_csv_by_first_column \$(rename_if_same_names \$first_file) $output_name

    other_files=\$(echo $files_to_join | cut -d' ' -f2-)
    for file in \$other_files; do
        sorted_file=\$(echo \$file |sed "s/.csv/_sorted.csv/g")   
        sort_csv_by_first_column \$file \$sorted_file
        join -t',' -1 1 -2 1 \$(rename_if_same_names $output_name) \$sorted_file > $output_name
    done
    """

    stub:
    """
    touch ${output_name}
    """
}
