process FILTER_OMARK_CONTAMINATIONS {
    tag "${annotation_method}/${id}"
    // label 'agat'
    scratch false
    memory '20 GB'
    stageInMode 'copy'
    
    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), path(fasta), path(ump), val(threshold)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations_c_${threshold}.txt"), emit: contaminations
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations_filtered_c_${threshold}.fa"), emit: fasta
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations_filtered_c_${threshold}.gff"), emit: gff

    script:
    def gff_for_omark = "${id}_${annotation_method}_for_omark.gff"
    def contaminations = "${id}_${annotation_method}_contaminations_c_${threshold}"
    def contaminations_filtered = "${id}_${annotation_method}_contaminations_filtered_c_${threshold}"
    """
    module unload nextflow
    ~/.conda/envs/python_plots/bin/python3 /shared/projects/metainvert/exogap/bin/reformat_gff_for_omark.py -i $gff -o $gff_for_omark

    # identify contamination and generation of proteins's fasta
    ~/.conda/envs/python_plots/bin/python3 /shared/projects/metainvert/exogap/bin/contamination_chromosome_filtering.py -i "." -g $gff_for_omark -f $fasta -o $contaminations -t $threshold
    mv ${contaminations}.fa ${contaminations_filtered}.fa

    # generation of filtered gff
    module load agat
        awk '
    /^-- Contaminant genes --[[:space:]]*\$/ {
        in_section=1
        skip_header=1
        next
    }

    /^-- .* --[[:space:]]*\$/ && in_section {
        in_section=0
        next
    }

    in_section {
        if (skip_header) { skip_header=0; next }
        if (\$1 != "") print \$1
    }
    ' ${contaminations}.txt | sort -u > contaminant_genes.txt

    agat_sp_filter_feature_from_keep_list.pl --gff $gff --keep_list contaminant_genes.txt --output ${contaminations_filtered}.gff
    """

    stub:
    """
    touch ${contaminations}.txt  ${contaminations}.fa  ${contaminations}_filtered.gff
    """
}
