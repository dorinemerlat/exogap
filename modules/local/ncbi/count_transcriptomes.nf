process COUNT_TRANSCRIPTOMES {
    tag "COUNT_TRANSCRIPTOMES_${name}"

    publishDir "out/data/transcriptomes/${taxid}/*", mode: 'copy'

    cache 'lenient'

    conda (params.enable_conda ? 'entrez-direct==16.2--he881be0_1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ncbi/edirect:12.5':
        'ncbi/edirect:12.5' }"

    input:
    tuple val(taxid), val(name), val(parents_count)

    output:
    tuple val(taxid), val(name), path{"${taxid}_transcriptomes.list"}, path{"${taxid}_transcriptomes_count.list"}, path{"${taxid}_transcriptomes_taxids.list"}

    script:
    """
    export NCBI_API_KEY=${params.ncbi_api_key}

    { # try a first time
        esearch -db nuccore -query '(txid${taxid}[Organism:exp]) AND "tsa master"[Properties]' > esearch.out

    } || { # try a second time
        sleep \$(shuf -i 5-30 -n 1)
        esearch -db nuccore -query '(txid${taxid}[Organism:exp]) AND "tsa master"[Properties]' > esearch.out
    }

    count=\$(grep "<Count>" esearch.out |cut -d '>' -f 2 |cut -d '<' -f1)

    if [[ \$count != 0 ]] ; then
        timer=\$((3 * (1 + \$RANDOM % $parents_count)))
        sleep \${timer}
        efetch -format gb < esearch.out > efetch.out
        grep "^ACCESSION" efetch.out |awk '{ print \$2}' > ${taxid}_transcriptomes.list
        grep -o '/db_xref="taxon:[0-9]*"' efetch.out | cut -d':' -f2 | cut -d'"' -f1 > ${taxid}_transcriptomes_taxids.list

        transcriptomes_count=\$(wc -l < ${taxid}_transcriptomes.list)
        echo "${taxid},${name},\${transcriptomes_count}" > ${taxid}_transcriptomes_count.list

    else
        touch ${taxid}_transcriptomes.list
        touch ${taxid}_transcriptomes_taxids.list
        echo "${taxid},${name},0" > ${taxid}_transcriptomes_count.list
    fi
    """

    // stub:
    // """
    //     touch ${taxid}_transcriptomes.list
    //     echo "${taxid},${name},0" > ${taxid}_transcriptomes_count.list
    // """
}
