process REPEATLANDSCAPE {
    tag "REPEATLANDSCAPE_${meta.id}"

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'quay.io/biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(align)

    output:
    tuple val(meta), path("${align}.html")

    script:
    """
    # Allow to find calcDivergenceFromAlign.pl and createRepeatLandscape.pl
    export PERL5LIB=/usr/local/share/RepeatMasker/

    # Re-calculate divergence adapted to GC content
    calcDivergenceFromAlign.pl -s ${align}.divsum -a ${align}.with_div $align

    # Repeat landscape
    createRepeatLandscape.pl -div ${align}.divsum -g ${meta.genome_size} > ${align}.html
    """
}
