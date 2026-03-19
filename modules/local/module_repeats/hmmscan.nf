// HMMSCAN: run rRNA HMM scan on the "for_mchelper" FASTA produced above
// This process is expected to consume the tuple emitted as `for_mchelper`
process HMMSCAN {
    tag "${id}"
    cpus 5
    label 'hmmer'
    time '2h'

    input:
    tuple val(id), val(molecule), path(fasta), path(hmm_db)

    output:
    // Save the raw hmmscan table and a filtered table with the grep selection
    tuple val(id), val(molecule), path("${id}_tes_vs_${molecule}.hmm")

    script:
    """
    hmmpress ${hmm_db}
    hmmscan --tblout ${tblout} -E 10 --noali --cpu ${task.cpus} "${id}_tes_vs_${molecule}.hmm" ${fasta}
    """
}
