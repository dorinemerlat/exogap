/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowExogap.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { fromSamplesheet } from 'plugin/nf-validation'

genomes = Channel.fromSamplesheet("input", skip_duplicate_check: true)
            .map { meta, fasta -> [ meta.name.toLowerCase().replace(' ', '-'), meta, file(fasta) ] }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// include { INPUT_CHECK                       } from '../subworkflows/input_check'
include { GET_INFORMATIONS_ABOUT_GENOMES    } from '../subworkflows/local/preprocess/get_informations_about_genomes'
include { PREPARE_GENOMES                   } from '../subworkflows/local/preprocess/prepare_genomes'
// include { ANALYSE_GENOME_QUALITY            } from '../subworkflows/preprocess/analyse_genome_quality'
// include { ANNOTATE_REPEATS                  } from '../subworkflows/annotate_repeats'
// include { ANNOTATE_PROTEIN_CODING_GENES     } from '../subworkflows/annotate_protein_coding_genes'
// include { ANNOTATE_NON_CODING_GENES         } from '../subworkflows/annotate_non_coding_genes'
// include { POSTPROCESS                       } from '../subworkflows/post_process'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def gather_genomes(genomes) {
    return genomes
        .map{it -> [it[0], it[1], it[2]]} // remove extra fiels
        .map { id, meta, file -> [id, ["${meta.taxid}": meta.name], file] } // keep only the IDs and essential meta data
        .flatten() // flat and list to get a destructed channel
        .toList()
        .map { genome -> genome.withIndex().collect { it, index -> [index % 3, it] } } // group by two: all the IDs, and all the meta data
        .flatMap { it } // have an item for IDs and one for meta data
        .groupTuple() // group to
        .map { index, it -> it.unique().sort() }
        .toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow EXOGAP {

    ch_versions = Channel.empty()

    // get taxonomy, datasets...
    GET_INFORMATIONS_ABOUT_GENOMES(genomes)

    // preprocess genomes
    PREPARE_GENOMES(GET_INFORMATIONS_ABOUT_GENOMES.out.genomes)

    // GET_INFORMATIONS_ABOUT_GENOMES(genomes)
    // PREPARE_GENOMES(genomes)
    // // ... .set { genomes } // add genomes size from PREPARE_GENOMES in the genomes of GET_INFORMATIONS_ABOUT_GENOMES
    // ANALYSE_GENOME_QUALITY(genomes)

    // //
    // // SUBWORKFLOW: ANNOTATE_REPEATS
    // //

    // // execute repeats annotation
    // if (params.annotate_repeats) {
    //     // ANNOTATE_REPEATS(genomes)
    //     //     REPETITIVE_ELEMENTS(PREPROCESS_GENOMES.out.genomes, PREPROCESS_GENOMES.out.newick)

    //    // -> masked_genomes
    // }
    // else {
    //     // ... .set{ genomes } // use masked genomes
    //     // integrate gff
    //     // ... .gff

    //     // -> masked_genomes
    // }

    // // ANNOTATE_PROTEIN_CODING_GENES(masked_genomes)
    // // ANNOTATE_NON_CODING_GENES(masked_genomes)

    // // POSTPROCESS(ANNOTATE_PROTEIN_CODING_GENES.out, ANNOTATE_NON_CODING_GENES.out)
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
