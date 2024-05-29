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
    .map { id, meta, fasta -> [ id, Utils.updateLinkedHashMap(meta, 'repeats_gff',   ['personal_set':meta.repeats_gff]), fasta ] }
    .set { genomes }

// bank blast
if (params.blastdb_local) {
    Channel.fromPath(params.blastdb_local + "/*")
        .collect()
        .map { it -> [it.find { it =~ ".phr" }.toString().replaceFirst(/.phr/, ""), it] }
        .set { blastdb }

} else if (params.blastdb_to_download ) {
    if (params.blastdb_to_download in ["nr", "swissprot", "refseq_prot"]) {
        Channel.from(params.blastdb).set { blastdb }
    } else {
        println "Invalid blastdb_to_download value. Please specify 'nr', 'swissprot', or 'refseq_prot'."
        exit 1
    }

} else if (params.blastdb_local && params.blastdb_to_download) {
    println "Please specify either 'blastdb_local' or 'blastdb_to_download', not both."
    exit 1

} else {
    println "Please specify either 'blastdb_local' or 'blastdb_to_download'."
    exit 1
}

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
include { ANALYSE_GENOME_QUALITY        } from '../subworkflows/local/analyse_genome_quality'
include { GET_DATASETS                  } from '../subworkflows/local/get_datasets'
include { GET_BLASTDB                   } from '../subworkflows/local/get_blastdb'
include { ANNOTATE_REPETITIVE_ELEMENTS  } from '../subworkflows/local/annotate_repetitive_elements'
include { ANNOTATE_PROTEIN_CODING_GENES } from '../subworkflows/local/annotate_protein_coding_genes'
include { ANNOTATE_NON_CODING_GENES     } from '../subworkflows/local/annotate_non_coding_genes'
include { POSTPROCESS                   } from '../subworkflows/local/postprocess'


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

    // analyse genome quality 
    ANALYSE_GENOME_QUALITY(PREPARE_GENOMES.out.genomes, PREPARE_GENOMES.out.newick, PREPARE_GENOMES.out.info)

    // // download datasets
    // GET_DATASETS(PREPARE_GENOMES.out.genomes)

    // if (params.blastdb_to_download) {
    //     GET_BLASTDB(blastdb)
    //     GET_BLASTDB.out.blastdb.set { blastdb }
    // }

    // // repetitive elements annotation
    // PREPARE_GENOMES.out.genomes
    //     .map { id, meta, fasta -> [ id, Utils.updateLinkedHashMap(meta, 'repeats_gff', meta.repeats_gff + ['repeatmasker': null]), fasta ] }
    //     .branch {
    //         no_repeats: it[1].repeats_gff.personal_set == null
    //         with_repeats: it[1].repeats_gff.personal_set != null }
    //     .set { genomes_for_repeats }

    // ANNOTATE_REPETITIVE_ELEMENTS(genomes_for_repeats.no_repeats)

    // // define the gff to use (anotate_repetitive_elements output or personal gff)
    // ANNOTATE_REPETITIVE_ELEMENTS.out.genomes
    //     .concat(genomes_for_repeats.with_repeats)
    //     .map { id, meta, fasta -> [ id, Utils.updateLinkedHashMap(meta, 'repeats_gff',  meta.repeats_gff + ['dataset': [meta.repeats_gff.personal_set, meta.repeats_gff.repeatmasker].findAll{ it != null }[0]]), fasta ] }
    //     .set { masked_genomes }

    // GET_DATASETS.out.genomes
    //     .join(masked_genomes)
    //     .map { id, meta1, fasta1, meta2, fasta2 -> [id, meta1 + ['repeats_gff': meta2.repeats_gff], fasta2]}
    //     .set { genomes_for_annotation }

    // ANNOTATE_PROTEIN_CODING_GENES( genomes_for_annotation, blastdb )

    // ANNOTATE_NON_CODING_GENES(genomes_for_annotation)

    // POSTPROCESS(genomes_for_annotation,
    //             ANNOTATE_PROTEIN_CODING_GENES.out.gff,
    //             ANNOTATE_NON_CODING_GENES.out.infernal,
    //             ANNOTATE_NON_CODING_GENES.out.barrnap_nucl,
    //             ANNOTATE_NON_CODING_GENES.out.barrnap_mito,
    //             ANNOTATE_NON_CODING_GENES.out.rnammer,
    //             genomes_for_annotation.map { id, meta, fasta -> [id, meta, meta.repeats_gff.dataset]})

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
