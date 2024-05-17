process BUSCO {
    tag "${id}_${dataset}_${mode}"
    cpus 30
    time '1d'
    label 'busco'
    stageInMode 'copy'
    stageOutMode 'copy'

    input:
    tuple val(id), val(meta), path(genome), val(dataset), val(mode)

    output:
    tuple val(id), path("busco_${id}_${dataset}_${mode}.json"),     emit: json
    tuple val(id), path("busco_${id}_${dataset}_${mode}.txt"),      emit: txt
    tuple val(dataset), path("busco_${id}_${dataset}_${mode}.csv"), emit: csv

    script:
    """
    busco -i $genome -l $dataset -o $id -m $mode -c $task.cpus
    mv ${id}/short_summary.*.json busco_${id}_${dataset}_${mode}.json
    mv ${id}/short_summary.*.txt busco_${id}_${dataset}_${mode}.txt

    jq -r '.results | [.Complete, ."Single copy", ."Multi copy", .Fragmented, .Missing, .n_markers] | @csv' busco_${id}_${dataset}_${mode}.json \
        | awk -v id=$id -v mode=$mode -v dataset=$dataset 'BEGIN {print "id,mode,dataset,complete,single_copy,multi_copyfragmented,missing,n_markers"} {print id "," mode "," dataset "," \$0}' \
        > busco_${id}_${dataset}_${mode}.csv
    """

    stub:
    """
    touch busco_${id}_${dataset}_${mode}.json
    touch busco_${id}_${dataset}_${mode}.txt
    touch busco_${id}_${dataset}_${mode}.csv
    """
}