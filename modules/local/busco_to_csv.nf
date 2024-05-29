process BUSCO_TO_CSV {
    tag "${id}_${dataset}_${mode}"
    label 'jq'
    scratch false

    input:
    tuple val(id), val(dataset), val(mode), path(json)

    output:
    tuple val(dataset), path("busco_${id}_${dataset}_${mode}.csv")

    script:
    """
    jq -r '.results | [.Complete, ."Single copy", ."Multi copy", .Fragmented, .Missing, .n_markers] | @csv' $json \
        | awk -v id=$id -v mode=$mode -v dataset=$dataset 'BEGIN {print "id,mode,dataset,complete,single_copy,multi_copyfragmented,missing,n_markers"} {print id "," mode "," dataset "," \$0}' \
        > busco_${id}_${dataset}_${mode}.csv
    """

    stub:
    """
    touch busco_${id}_${dataset}_${mode}.csv
    """
}