process CLASSIFY_GENES_VS_REPEATS {
    tag "${annotation_method}/${id}"
    // label 'busco'
    scratch false
    stageInMode 'copy'
    conda 'envs/python_plots.yml'

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins), path(repeats)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_genes_vs_merged_repeats.tsv"), emit: merge
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_genes_vs_no_merged_repeats.tsv"), emit: no_merge

    script:
    """
    # show which python/conda are visible (debug)
    module unload nextflow
    module load bedtools
    command -v conda || true
    which python || true

    mkdir tmp

    conda run -n python_plots --no-capture-output python /shared/projects/metainvert/exogap/bin/classify_genes_vs_repeats.py \
        --proteins ${proteins} \
        --repeats ${repeats} \
        --output ${id}_${annotation_method}_genes_vs_merged_repeats.tsv \
        --merge-repeats \
        --tmp tmp

    conda run -n python_plots --no-capture-output python /shared/projects/metainvert/exogap/bin/classify_genes_vs_repeats.py \
        --proteins ${proteins} \
        --repeats ${repeats} \
        --output ${id}_${annotation_method}_genes_vs_no_merged_repeats.tsv \
        --tmp tmp
    """

    stub:
    """
    touch ${id}_${annotation_method}_genes_vs_merged_repeats.tsv ${id}_${annotation_method}_genes_vs_no_merged_repeats.tsv 
    """
} 