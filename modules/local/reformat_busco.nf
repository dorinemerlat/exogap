process REFORMAT_BUSCO {
    tag "REFORMAT_BUSCO_${id}_${dataset}_${mode}"

    input:
    tuple val (id), val(meta), path(genome), path(json), val(dataset), val(mode)

    output:
    tuple val (id), val(meta), path("busco_${id}.${dataset}.tsv")

    script:
    """
    jq -r ' [.C,.S,.D,.F,.M,.dataset_total_buscos] | @csv' $json \
    | awk -v header="genome_id,genome_name,taxid,dataset,mode,complete,single,duplicate,fragmented,missing,total" \
        -v dataset=$dataset -v id=$id -v name="${meta.name}" -v taxid=$meta.taxid -v mode=$mode -v OFS=',' \
        'BEGIN{OFS=","; print header} {print id, name, taxid, dataset, \$0}' \
    | sed "s/\\"//g" \
    > busco_${id}.${dataset}.tsv
    """

    stub:
    """
    touch "busco_${id}.${dataset}.tsv"
    """
}
