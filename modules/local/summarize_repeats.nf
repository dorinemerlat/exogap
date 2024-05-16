process SUMMARIZE_REPEATS {
    tag "SUMMARIZE_REPEATS"
    input:
    tuple val(meta), path(out)


    output:
    tuple val(meta), path("${meta.id}.csv")

    script:
    """
	sed "s|SINE?|SINE|g; s|RC/|DNA/|g; s|Simple_repeat|Satellite/Simple_repeat|g; s|ARTEFACT|Other/Artefact|g; s|Integrated-virus|Other/Integrated_virus|g" $out > temp.csv

	awk -F',' 'BEGIN {OFS=FS} NR==1 {print \$1, \$2, \$3, "te_class", "te_family", "te_length_total"} NR>1 {te_length[\$5]+=\$8} END {for (te_name in te_length) {split(te_name,te_name_parts,"/"); count[te_name]++; print \$1, \$2, \$3, te_name_parts[1], te_name_parts[2], te_length[te_name], count}}' $temp.csv > ${meta.id}.csv

    plot_repeats.R ${meta.id}.csv
    """

    stub:
    """
    touch ${meta.id}.csv
    """
}
