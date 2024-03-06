process DOWNLOAD_TSA {
    tag "DOWNLOAD_TSA_${clade_name}"

    publishDir "out/data/sra/${clade_name}/*"

    cache 'lenient'

    conda (params.enable_conda ? 'sra-tools==3.0.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.0.8--h9f5acd7_0' :
        'biocontainers/sra-tools:3.0.8--h9f5acd7_0' }"

    input:
    tuple val(clade_taxid), val(clade_name), path(id_list)

    output:
    tuple val(clade_taxid), val(clade_name), path("${clade_taxid}_transcriptomes.fasta")

    script:
    """
    prefetch --option-file $id_list

    srr_directories=\$(find . -type f ! -name '.*')
    for srr in \${srr_directories} ; do
        fasterq-dump --fasta \${srr}
    done

    cat *fasta > ${clade_taxid}_transcriptomes.fasta
    """

    // stub:
    // """
    //     touch ${taxid}_transcriptomes.list
    //     echo "${taxid},${name},0" > ${taxid}_transcriptomes_count.list
    // """
}
