process EXTRACT_SEQUENCES {
    tag "EXTRACT_SEQUENCES_${meta.id}"
    publishDir "${params.outdir}/results/${meta.id}/repetitive-elements"

    conda (params.enable_conda ? 'agat==1.2.0--pl5321hdfd78af_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.0--pl5321hdfd78af_0':
        'quay.io/biocontainers/agat:1.2.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(genome)
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${genome}-repetitive-elements.fa")

    script:
    """
    agat_sp_extract_sequences.pl -g $gff -f $genome -t dispersed_repeat -o
    """
}
