// HMMSCAN: run rRNA HMM scan on the "for_mchelper" FASTA produced above
// This process is expected to consume the tuple emitted as `for_mchelper`
process HMMSCAN {
    tag "${id}"
    cpus 5
    time '2h'

    input:
    tuple val(id), val(meta), path(fasta), path(hmm_db), val(tblout)

    output:
    // Save the raw hmmscan table and a filtered table with the grep selection
    tuple val(id), val(meta), path(tblout)

    script:
    """
    module load hmmer
    hmmpress ${hmm_db}
    hmmscan --tblout ${tblout} -E 10 --noali --cpu ${task.cpus} ${hmm_db} ${fasta}
    """
}
