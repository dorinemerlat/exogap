process GFF_TO_TSV {
    tag "$id"
    // label 'exogap_python'
    memory { (5.GB * (task.attempt * task.attempt)) }
    maxRetries 5

    input:
    tuple val(id), val(meta), path(gff)

    output:
    tuple val(id), val(meta), path("${id}_repeats_summary.tsv")

    script:
    """
    /home/merlat/.conda/envs/biopython/bin/python /tempor/merlat/exogaptwo/bin/gff_to_tsv.py -i $gff -o ${id}_repeats_summary.tsv
    """

    stub:
    """
    touch ${id}_repeats_summary.tsv
    """
}
