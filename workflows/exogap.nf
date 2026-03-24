/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPROCESSING } from '../subworkflows/local/preprocessing.nf'
include { REPEATS_ANNOTATION } from '../subworkflows/local/repeats_annotation.nf'
include { GENES_ANNOTATION } from '../subworkflows/local/genes_annotation.nf'
include { NCGENES_ANNOTATION } from '../subworkflows/local/ncgenes_annotation.nf'
include { POSTPROCESSING } from '../subworkflows/local/postprocessing.nf'

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
            REPEATS_ANNOTATION(PREPROCESSING.out.genomes, PREPROCESSING.out.newick, PREPROCESSING.out.genome_stats)

        }

        // if (params.module_genes == true) {
        //     print "Running module gene annotation"

        //     if (params.module_repeats == true) {
        //         REPEATS_ANNOTATION.out.unmasked_genomes.set { unmasked_genomes }
        //         REPEATS_ANNOTATION.out.masked_genomes.set { masked_genomes }
        //         REPEATS_ANNOTATION.out.repeat_gff.set { repeat_gff }
        //     } else {
        //         PREPROCESSING.out.genomes.set { unmasked_genomes }
        //         PREPROCESSING.out.masked_genomes.set { masked_genomes }
        //         Channel.empty().set { repeat_gff }
        //     }

        //     GENES_ANNOTATION( unmasked_genomes, masked_genomes, repeat_gff )
        // }

        // if (params.module_ncgenes == true) {
        //     print "Running module ncRNA annotation"

        //     if (params.module_repeats == true) {
        //         REPEATS_ANNOTATION.out.unmasked_genomes.set { genomes_for_ncgenes_annotation }
        //     } else {
        //         PREPROCESSING.out.genomes.set { genomes_for_ncgenes_annotation }
        //     }

        //     NCGENES_ANNOTATION( genomes_for_ncgenes_annotation )
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
