process GFF_TO_TSV {
    tag "$id"
    label 'exogap_python'
    memory { (5.GB * (task.attempt * task.attempt)) }
    maxRetries 5
    
    input:
    tuple val(id), val(meta), path(gff)

    output:
    tuple val(id), val(meta), path("${id}_repeats.tsv")

    script:
    """
    gff_to_tsv.py -i $gff -o ${id}_repeats.tsv.tmp
    awk -v id="$id" 'BEGIN {OFS="\\t"} NR == 1 {print "genome", \$0} NR > 1 {print id, \$0}' ${id}_repeats.tsv.tmp > ${id}_repeats.tsv
    """

    stub:
    """
    touch ${id}_repeats.tsv
    """
}
