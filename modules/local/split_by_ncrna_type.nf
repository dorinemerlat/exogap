process SPLIT_BY_NCRNA_TYPE {
    scratch true
    tag "SPLIT_BY_NCRNA_TYPE_${id}"

    input:
    tuple val(id), val(meta), path(genome), path(gff), path(families)

    output:
    tuple val(id), val(meta), path(genome), path("${id}_*.gff", arity: '1..*')

    script:
    """
    while IFS= read -r line; do
        if [[ "\${line:0:1}" != "#" ]] && [[ "\$line" =~ ";mdlaccn=" ]]; then
            # Extract the pattern from the line
            pattern=\$(echo "\$line" | awk -F';mdlaccn=' '{print \$2}' | awk -F';' '{print \$1}')

            # Search for the pattern in rfam_class file and extract the family
            family=\$(grep "\$pattern" $families | awk '{print \$2}')

            # Append family information to the line and save to output file
            echo "\$line;family=\$family" >> ${id}_families.gff
        fi
    done < $gff

    for i in \$(grep -oP '(?<=;family=)[^;]*\$' ${id}_families.gff |sort -u); do
        grep -o 'family='\$i > ${id}_\${i}.gff
    done
    """

    stub:
    """
    touch ${id}_rRNA.gff ${id}_snRNA.gff
    """
}
