process GENERATE_RUNLIST {
    tag "${id}"
    scratch false

    input:
    tuple val(id), val(meta)

    output:
    tuple val(id), val(meta), path("Runlist_${id}.txt")

    script:
    def sra_list = meta.SRA ? meta.SRA.replaceAll(';', ' ') : null

    if (sra_list) {
        """
        generate_runlist.sh ${sra_list} > Runlist_${id}.txt
        """
    } else {
        """
        # No SRA provided
        touch Runlist_${id}.txt
        """
    }

    stub:
    """
    touch Runlist_${id}.txt
    """
}
