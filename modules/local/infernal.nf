process INFERNAL {
    tag "INFERNAL_${id}"
    label 'infernal'
    cpus 40

    input:
    tuple val(id), val(meta), path(genome), path(rfam_clanin), path(rfam_cm), path(family)

    output:
    tuple val(id), val(meta), path("${id}_infernal.gff")

    script:
    """
    residues=\$(esl-seqstat $genome | grep 'Total # residues' ${id}_infernal.stats | cut -d ':' -f 2)
    Z=\$(echo "\$residues*2/1000000" |bc -l)

    cmscan -Z \$Z --cpu $task.cpus --cut_ga --rfam --nohmmonly --tblout ${id}.tblout -o ${id}.out \
            --fmt 2 --clanin $rfam_clanin $rfam_cm $genome

    infernal-tblout2gff.pl --cmscan --fmt2 --all ${id}.tblout > ${id}_infernal.gff.tmp

    grep -v '^#' ${id}_infernal.gff.tmp \
        | awk '{ nb++; print \$1"\t"molecule_type\tgene\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\tID=infernal_"nb";Name="\$2";"\$9"\n" \
            \$1"\tmolecule_type\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\tID=infernal_"nb"-rRNA;Parent=infernal_"nb";Name="\$2";"\$9}' \
        > ${id}_infernal.gff.tmp.tmp

    while IFS= read -r line; do
        if [[ "\${line:0:1}" != "#" ]] && [[ "\$line" =~ ";mdlaccn=" ]]; then
            # Extract the pattern from the line
            pattern=\$(echo "\$line" | awk -F';mdlaccn=' '{print \$2}' | awk -F';' '{print \$1}')

            # Search for the pattern in rfam_class file and extract the family
            family=\$(grep "\$pattern" $family | awk '{print \$2}')
            molecule_type=\$(echo \$family |awk -F',' '{if (\$1 == "Cis-reg") print \$1; else if (\$3 == "snoRNA") print \$3; else if (\$3 == "splicing") print \$2; else print \$NF}' )

            # Append family information to the line and save to output file
            echo "\$line;family=\$family" | sed "s/molecule_type/\$molecule_type/g" >>  ${id}_infernal.gff
        fi
    done < ${id}_infernal.gff.tmp.tmp
    """

    stub:
    """
    touch ${id}_infernal.gff
    """
}
