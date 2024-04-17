process UPDATE_LOCALDB {
    tag "UPDATE_LOCALDB_${id}"
    cpus 10

    input:
    val(blastdb)

    output:
    tuple path("${blastdb}*")

    script:
    """
    update_blastdb.pl --num_threads $task.cpus --source ncbi --decompress ${blastdb}
    """

    stub:
    """
    for i in "pal", "pdb", "phr", "pin", "pjs", "pot", "psq", "ptf", "pto" do
        touch ${id}_db.\${i}
    done
    """
}

