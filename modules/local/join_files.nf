process JOIN_FILES {
    scratch true
    // in tag, replace characters: [, ] by nothing
    tag "GATHER_FILES_${meta.name}"

    input:
    tuple path(file, stageAs: "file1.csv"), path(new_file), val(output_name)

    output:
    path output_name

    script:
    """
    sort -t',' -k1,1 file1.csv > sorted_file1.csv
    sort -t',' -k1,1 $new_file > sorted_file2.csv

    join -t',' -1 1 -2 1 sorted_file1.csv sorted_file2.csv > $output_name
    """

    stub:
    """
    touch ${output_name}
    """
}