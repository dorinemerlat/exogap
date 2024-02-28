
class Tools {
    static gather_genomes(List a) {
        println(a)
    }
    // static List gather_genomes(List genomes) {
    //     return genomes
    //         .map{it -> [it[0], it[1], it[2]]} // remove extra fiels
    //         .map { id, meta, file -> [id, ["${meta.taxid}": meta.name], file] } // keep only the IDs and essential meta data
    //         .flatten() // flat and list to get a destructed channel
    //         .toList()
    //         .map { genome -> genome.withIndex().collect { it, index -> [index % 3, it] } } // group by two: all the IDs, and all the meta data
    //         .flatMap { it } // have an item for IDs and one for meta data
    //         .groupTuple() // group to
    //         .map { index, it -> it.unique().sort() }
    //         .toList()
    // }
}
