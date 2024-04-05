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
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { fromSamplesheet } from 'plugin/nf-validation'

Channel.fromSamplesheet("input", skip_duplicate_check: true)
    .map { meta, fasta -> [ meta.name.toLowerCase().replace(' ', '-'), meta, file(fasta) ] }
    .map { id, meta, fasta -> [ id, Utils.updateLinkedHashMap(meta, 'main_protein_set', ['personal_set': meta.main_protein_set]), fasta ] }
    .map { id, meta, fasta -> [ id, Utils.updateLinkedHashMap(meta, 'training_protein_set', ['personal_set': meta.training_protein_set]), fasta ] }
    .map { id, meta, fasta -> [ id, Utils.updateLinkedHashMap(meta, 'transcript_set', ['personal_set': meta.transcript_set]), fasta ] }
    .map { id, meta, fasta -> [ id, Utils.updateLinkedHashMap(meta, 'repeats_gff',  meta.repeats_gff), fasta ] }
    .set { genomes }

// bank blast
Channel.fromPath(params.blast_db_local + "/*")
    .collect()
    .map { it -> [it.find { it =~ ".phr" }.toString().replaceFirst(/.phr/, ""), it] }
    .set { blast_db }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// include { INPUT_CHECK                       } from '../subworkflows/input_check'
include { PREPARE_GENOMES               } from '../subworkflows/local/prepare_genomes'
// include { ANALYSE_GENOME_QUALITY            } from '../subworkflows/preprocess/analyse_genome_quality'
include { GET_DATASETS                  } from '../subworkflows/local/get_datasets'
include { ANNOTATE_REPETITIVE_ELEMENTS  } from '../subworkflows/local/annotate_repetitive_elements'
include { ANNOTATE_PROTEIN_CODING_GENES } from '../subworkflows/local/annotate_protein_coding_genes'
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

    // get taxonomy and preprocess genomes
    PREPARE_GENOMES(genomes)

    genomes.branch {
        no_repeats: it[1].repeats_gff == null
        with_repeats: other }
        .set { genomes }

    if (genomes.no_repeats.count() != 0) {
        ANNOTATE_REPETITIVE_ELEMENTS(genomes.no_repeats)

        ANNOTATE_REPETITIVE_ELEMENTS.out.genomes
            .concat ( genomes.with_repeats )
            .set { genomes }
    }

    // execute genome annotation
    GET_DATASETS(PREPARE_GENOMES.out.genomes)

    GET_DATASETS.out.genomes.join(genomes)
        .map { id, meta1, fasta1, meta2, fasta2 -> [id, Utils.updateLinkedHashMap(meta1, 'repeats_gff', meta2.repeats_gff), fasta1]}
        .set { genomes }

    ANNOTATE_PROTEIN_CODING_GENES( genomes, blast_db )
    // // ANNOTATE_NON_CODING_GENES(masked_genomes)

    // // POSTPROCESS(ANNOTATE_PROTEIN_CODING_GENES.out, ANNOTATE_NON_CODING_GENES.out)

    // emit:
    //     repetitive_elements = ANNOTATE_REPETITIVE_ELEMENTS.out.masked
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
