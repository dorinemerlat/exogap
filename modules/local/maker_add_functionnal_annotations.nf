process MAKER_ADD_FUNCTIONNAL_ANNOTATIONS {
    tag "MAKER_ADD_FUNCTIONNAL_ANNOTATIONS_${id}_${iteration}"
    label 'maker'

    input:
    tuple val(id), val(meta), path(genome), path(train), path(test), path(specie), val(iteration)

    output:
    tuple val(id), val(meta), path(genome), path(train), path(test), path("config/species/${id}_*"), val(iteration)

    script:
    """
    maker_functional_gff uniprot_sprot.db $blastp $gff > hsap_contig.renamed.putative_function.gff
    maker_functional_fasta uniprot_sprot.db $blastp $proteins > hsap_contig.maker.proteins.renamed.putative_function.fasta
    maker_functional_fasta uniprot_sprot.db $blastp $transcripts > hsap_contig.maker.transcripts.renamed.putative_function.fasta
    """

    stub:
    """
    mkdir -p config/species/
    for i in $specie; do
        cp \$i config/species/
    done
    """
}
