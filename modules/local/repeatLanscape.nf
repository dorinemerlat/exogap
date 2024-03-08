process repeatlandscape {
    tag "repeatlandscape_${meta.id}"
    label 'repeatmasker'

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
