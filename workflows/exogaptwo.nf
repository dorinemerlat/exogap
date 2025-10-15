/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPROCESSING } from '../subworkflows/local/module_preprocessing/preprocessing.nf'
include { REPEATS_ANNOTATION } from '../subworkflows/local/module_repeats/repeats_annotation.nf'
include { GENES_ANNOTATION } from '../subworkflows/local/module_genes/genes_annotation.nf'
include { NCGENES_ANNOTATION } from '../subworkflows/local/module_ncgenes/ncgenes_annotation.nf'
include { POSTPROCESSING } from '../subworkflows/local/module_postprocessing/postprocessing.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EXOGAPTWO {
    take:
        genomes

    main:
        print "Running workflow EXOGAP"

        print "Running module cleaning genomes"
        PREPROCESSING(genomes)

        if (params.module_repeats == true) {
            print "Running module repeat annotation"
            REPEATS_ANNOTATION(PREPROCESSING.out.genomes)
        }

        // if (params.module_genes == true) {
        //     print "Running module gene annotation"
        //     GENES_ANNOTATION()
        // }

        // if (params.module_ncgenes == true) {
        //     print "Running module ncRNA annotation"
        //     NCGENES_ANNOTATION()
        // }

        // print "Running module postprocessing"
        // POSTPROCESSING()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
