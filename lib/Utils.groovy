//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowBroadcast
import nextflow.script.ChannelOut

class Utils {

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        // This channel list is ordered by required channel priority.
        def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
        def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

        // Check that they are in the right order
        def channel_priority_violation = false
        def n = required_channels_in_order.size()
        for (int i = 0; i < n - 1; i++) {
            channel_priority_violation |= !(channels.indexOf(required_channels_in_order[i]) < channels.indexOf(required_channels_in_order[i+1]))
        }

        if (channels_missing | channel_priority_violation) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/\n" +
                "  The observed channel order is \n" +
                "  ${channels}\n" +
                "  but the following channel order is required:\n" +
                "  ${required_channels_in_order}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }

    //
    // Function to update a LinkedHashMap with a given key-value pair
    //
    public static LinkedHashMap updateLinkedHashMap(LinkedHashMap dict, Object key, Object value) {
        dict[key] = value
        return dict
    }

    public static String getBaseNameOrNull(String value) {
        if (value == null) {
            return null
        } else {
            return value.getBaseName()
        }
    }

    public static DataflowBroadcast gatherGenomes(DataflowBroadcast genomes) {
        return genomes
            .map { it -> [it[0], it[1], it[2]] } // remove extra fields
            .map { id, meta, file -> [id, ["${meta.taxid}": meta.name], file] } // keep only the IDs and essential meta data
            .flatten() // flatten the list to get a destructed channel
            .toList()
            .map { genome -> genome.withIndex().collect { it, index -> [index % 3, it] } } // group by two: all the IDs, and all the meta data
            .flatMap { it } // have an it for IDs and one for meta data
            .groupTuple() // group tuples
            .map { index, it -> it.unique().sort() }
            .toList()
    }

    public static DataflowBroadcast indexGenomesByLineage(DataflowBroadcast genomes) {
        return genomes
            .map { id, meta, fasta -> [meta.lineage.taxid, id, meta, fasta ] }
            .transpose()
    }


    public static List createEmptySet() {
        return [[null, ['taxid': null, 'count': null, 'other': null], null]]
    }
}

