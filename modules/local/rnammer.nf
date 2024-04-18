process RNAMMER {
    scratch true
    tag "RNAMMER_${id}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}_rnammer.gff")

    script:
    """
    rnammer -S euk -m tsu,ssu,lsu -multi -gff ${id}_rnammer.gff.tmp -f ${id}_rnammer.fa.tmp -h ${id}_rnammer.report < $genome

    # correct gff file
    grep -v '^#' ${id}_rnammer.gff.tmp \
        | awk -F'\t' -v OFS='\t' '\$9 == "8s_rRNA" { \$9 = "5s_rRNA" } { print }' \
        | awk '{ nb++; print \$1"\t"\$2"\t""gene""\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\tID=rnammer_"nb";Name="\$9";family=Gene,rRNA\n" \
            \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\tID=rnammer_"nb"-rRNA;Name="\$9";Parent=rnammer_"nb";family=Gene,rRNA"}' \
        > ${id}_rnammer.gff

    # correct fasta file
    sed "s|molecule=8s_rRNA|molecule=5s_rRNA|g" ${id}_rnammer.fa.tmp > ${id}_rnammer.fa

    # correct file
    sed -i 's/s_rRNA/S_rRNA/g' ${id}_rnammer.*
    """

    stub:
    """
    touch ${id}_rnammer.gff ${id}_rnammer.fa
    """
}
