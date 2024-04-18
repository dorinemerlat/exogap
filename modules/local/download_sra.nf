process DOWNLOAD_SRA {
    scratch true
    tag "DOWNLOAD_SRA_${specie_name}"
    cache 'lenient'
    label 'sratools'

    input:
    tuple val(clade_taxid), val(specie_name), val(specie_taxid), path(efetch_out), val(efetch_count)

    output:
    tuple val(clade_taxid), val(specie_name), val(specie_taxid), path("*_1.fastq"), path("*_2.fastq")

    script:
    """
    prefetch --option-file $efetch_out

    srr_directories=\$(find . -type d -name 'SRR*')

    for srr in \${srr_directories}; do
        fasterq-dump --split-files \${srr}
    done
    """

    stub:
    """
    touch sra_1.fastq sra_2.fastq
    """
}
