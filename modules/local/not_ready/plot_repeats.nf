process PLOT_REPEATS {
    input:
    tuple val(meta), val(genomes)

    output:
    path("repeat_overview.pdf")

    script:
    """
    plot_repeats.R
    """
}
