#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/exogaptwo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/exogaptwo
    Website: https://nf-co.re/exogaptwo
    Slack  : https://nfcore.slack.com/channels/exogaptwo
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// // Print help message if needed
// if (params.help) {
//     def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
//     def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
//     def String command = "nextflow run ${workflow.manifest.name} -c nextflow.config -profile singularity"
//     log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
//     System.exit(0)
// }

// // Validate input parameters
// if (params.validate_params) {
//     validateParameters()
// }
// print 'prout0'
// WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check is the user provided a samplesheet
if (params.input) {
    // Check if the samplesheet exists
    def samplesheet_file = file(params.input, checkIfExists: true)
    // Convert samplesheet to list of maps
    samplesheet = samplesheetToList(samplesheet_file, "$projectDir/assets/schema_input.json")
    log.info "Found ${samplesheet.size()} genomes in the samplesheet: ${params.input}"
} else {
    log.error "No samplesheet provided. You must provide a samplesheet with --input or in nextflow.config file to run the workflow."
    System.exit(1)
}

// check if at least one module is activated
if (params.module_repeats == false && params.module_genes == false && params.module_ncgenes == false) {
    log.error "You have deactivated all annotation modules. Please activate at least one module to run the workflow."
    System.exit(1)
}

// if module_repeats is activated, check its parameters
if (params.module_repeats == true) {
    // check download_famdb
    if (params.download_famdb != true && params.download_famdb != false) {
        log.error "Parameter error: --download_famdb must be true or false."
        System.exit(1)
    }
    // check famdb_path and famdb_to_download
    // if download_famdb is false, check that famdb_path is provided and if it is an existing directory containing dfam39_full.0.h5 and at least one other .h5 file with 1 to 16 in its name
    if (params.download_famdb == false) {
        if (params.famdb_path == null) {
            log.error "Parameter error: You set --download_famdb to false, so you must provide a local path to a FamDB library with --famdb_path."
            System.exit(1)
        } else {
            def famdb_dir = file(params.famdb_path, checkIfExists: true)
            if (!famdb_dir.isDirectory()) {
                log.error "Parameter error: --famdb_path must be a path to a directory containing a FamDB .h5 file."
                System.exit(1)
            } else {
                params.famdb_partitions = famdb_dir.listFiles().findAll { it.name ==~ /^dfam39_full(?:\.\d+)?\.h5$/ || it.name ==~ /^dfam39_(?:[1-9]|1[0-6])(?:\.\d+)?\.h5$/ }
                if (params.famdb_partitions.size() < 2) {
                    log.error "Parameter error: The directory provided in --famdb_path must contain at least the full FamDB (dfam39_full.0.h5) and one other partitioned FamDB (.h5 file with 1 to 16 in its name)."
                    System.exit(1)
                }
                log.info "Found ${params.famdb_partitions.size()} FamDB partition(s): ${params.famdb_partitions*.name.sort().join(', ')}"
                // do a list of the available partitions with the path of directory join to files names
            }
        }
    } else {
        // if download_famdb is true, check that famdb_to_download is valid
        if (params.famdb_to_download == null) {
            log.error "Parameter error: You set --download_famdb to true, so you must provide a FamDB ID to download from dfam.org with --famdb_to_download."
            System.exit(1)
        } else {
            if (!(params.famdb_to_download ==~ /^(?:[1-9]|1[0-6])$/)) {
                log.error "Parameter error: --famdb_to_download must be a string '1'..'16' corresponding to a FamDB library ID on dfam.org."
                System.exit(1)
            }
        }
    }
    // check group_consensus_sequences
    if (params.group_consensus_sequences != true && params.group_consensus_sequences != false) {
        log.error "Parameter error: --group_consensus_sequences must be true or false."
        System.exit(1)
    }

    // Vérification du paramètre optionnel reference_library
    if (params.reference_library) {
        def ref_lib_file = file(params.reference_library)
        if (!ref_lib_file.exists()) {
            log.warn "The file provided for --reference_library does not exist: ${params.reference_library}"
        } else {
            params.reference_library = ref_lib_file
            log.info "Reference library file provided for repeat annotation: ${params.reference_library}"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXTRACT GENOME INFORMATION FROM SAMPLESHEET
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

Channel.fromList(samplesheet)
    // create a unique ID for each genome based on its name, in lowercase and with spaces replaced by '-'
    .map { meta, fasta -> [ meta.name.toLowerCase().replace(' ', '-'), meta, file(fasta) ] }
    // make sure taxid is a string
    .map { id, meta, fasta -> [ id, Utils.updateTopLevelKey(meta, 'taxid', meta.taxid.toString()), fasta ] }
    // look if user has provided dataset to use
    .map { id, meta, fasta ->
        def updatedMeta = meta.collectEntries { key, value ->
            // search the meta parameter relative to the user dataset
            if (key.contains("_set")) {
                // look if user has provided a set for each specie...
                if (!value.isEmpty()) {
                    // ... if so, use it
                    [key, ["user_set": value]]
                } else {
                    // ... if not, look if a set is provided in the params for all species...
                    if (params.containsKey(key)) {
                        // ... if so, use it
                        [key, ["user_set": params."$key"]]
                    } else {
                        // ... if not, set it to null
                        [key, ["user_set": null]]
                    }
                // if user provides a set for each specie, use it
                }
            // for all other parameters, keep them as they are
            } else {
                [key, value]
            }
        }
        [id, updatedMeta, file(fasta)]
    }
    .set { genomes }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXOGAPTWO } from './workflows/exogaptwo'


// WORKFLOW: Run main nf-core/exogaptwo analysis pipeline

workflow {
    print "Starting workflow EXOGAP"
    EXOGAPTWO (genomes)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
