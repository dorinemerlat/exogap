process RENAME_REPEAT_FAMILIES {
    tag "$id"
    label 'bioawk'

    input:
    tuple  val(id), val(meta), path(families)

    output:
    tuple  val(id), val(meta), path("${id}_repeat_families.fa")

    script:
    def genome_code = meta.mnemonic.exogap
    """
    # Update header lines
    bioawk -v id="TE_${genome_code}" -c fastx '{split(\$name, arr, "#"); print ">"id"_"NR"#"arr[2]; print \$seq}' $families \
        | fold -w 60 \
        | sed -e 's/#unknown/#Unknown/' -e 's@/unknown@/Unknown@' \
        > ${id}_repeat_families.fa
    """

    stub:
    """
    touch ${id}_repeat_families.fa
    """
}
