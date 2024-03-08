process DOWNLOAD_SRA {
    tag "DOWNLOAD_SRA_${specie_name}"

    publishDir "out/data/sra/${specie_taxid}/*"

    cache 'lenient'

    conda (params.enable_conda ? 'sra-tools==3.0.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.0.8--h9f5acd7_0' :
        'biocontainers/sra-tools:3.0.8--h9f5acd7_0' }"

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

    // stub:
    // """
    //     touch ${taxid}_transcriptomes.list
    //     echo "${taxid},${name},0" > ${taxid}_transcriptomes_count.list
    // """
}
