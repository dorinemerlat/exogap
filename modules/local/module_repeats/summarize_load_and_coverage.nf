process SUMMARIZE_LOAD_AND_COVERAGE {
    tag "${id}"
    label 'exogap_tools'

    input:
    tuple val(id), val(meta), path(tsv)

    output:
    tuple val(id), val(meta), path("${id}_TE_load_and_coverage.csv")

    script:
    def bed_file = "temp.bed"
    def merged_bed = "merged.bed"
    """
    # Remove header, filter lines and extract columns 2, 5, and 6 (seqid, start, end)
    column_to_filter=\$(head ${tsv} | tr '\t' '\n' | nl | grep -P '\trepeat_type')

    tail -n +2 "${tsv}" \
        | awk -F'\\t' -v col=\$(head -n 1 "${tsv}" | tr '\\t' '\\n' | nl | grep -P '\\trepeat_type' | awk '{print \$1}') '{
            if (\$col == "Tandem Repeat" || \$col == "Transposable Element" || \$col == "Unknown") {
                print \$2, \$5, \$6
            }
        }' OFS='\\t' > "$bed_file"

    # Apply bedtools merge
    bedtools merge -i "$bed_file" > "$merged_bed"

    # Calculate total merged coverage
    total_length=\$(awk '{sum += (\$3 - \$2)} END {print sum}' "$merged_bed")
    coverage=\$(( 100 * \${total_length} / ${meta.assembly_size} ))

    # Count number of repeat elements (excluding header)
    load=\$((\$(wc -l < "${tsv}") - 1))

    # Write results to CSV
    echo "genome,coverage,load" > "${id}_TE_load_and_coverage.csv"
    echo "${id},\${coverage},\${load}" >> "${id}_TE_load_and_coverage.csv"
    """

    stub:
    """
    echo "genome,coverage,load" > "${id}_TE_load_and_coverage.csv"
    """
}