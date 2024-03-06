process SEARCH_SRA {
    tag "SEARCH_SRA_${specie_name}"

    cache 'lenient'

    conda (params.enable_conda ? 'entrez-direct==16.2--he881be0_1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ncbi/edirect:12.5':
        'ncbi/edirect:12.5' }"

    input:
    tuple val(clade_taxid), val(specie_name), val(specie_taxid), val(specie_count)

    output:
    tuple val(clade_taxid), val(specie_name), val(specie_taxid), path("efetch.out")

    script:
    """
    export NCBI_API_KEY=${params.ncbi_api_key}

    { # try a first time
         esearch -db sra -query '((((txid${specie_taxid}[Organism:exp]) AND "paired"[Layout]) AND "illumina"[Platform]) AND "rna data"[Filter]) AND "filetype fastq"[Properties]' \
            > esearch.out

    } || { # try a second time
        sleep \$(shuf -i 5-30 -n 1)
        esearch -db sra -query '((((txid${specie_taxid}[Organism:exp]) AND "paired"[Layout]) AND "illumina"[Platform]) AND "rna data"[Filter]) AND "filetype fastq"[Properties]' \
            > esearch.out
    }

    # remove header
    count=\$(grep "<Count>" esearch.out |cut -d '>' -f 2 |cut -d '<' -f1)

    if [[ \$count != 0 ]] ; then
        timer=\$((3 * (1 + \$RANDOM % $specie_count)))
        sleep \${timer}
        efetch -format runinfo < esearch.out |cut -f 1 -d ',' |sed 1d > efetch.out

    else
        touch efetch.out
    fi

    """

    // stub:
    // """
    //     touch ${taxid}_transcriptomes.list
    //     echo "${taxid},${name},0" > ${taxid}_transcriptomes_count.list
    // """
}
