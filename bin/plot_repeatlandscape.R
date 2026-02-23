#!/usr/bin/env Rscript

## ------------------------------------------------------------
## Repeat landscapes: stacked histograms per genome (by TE order)
## Output: one PDF/SVG containing all genomes in phylogeny order
## ------------------------------------------------------------

## (Optional) remove if Nextflow sets it
setwd("/tempor/merlat/exogaptwo")

## ------------------------------------------------------------
## Packages (minimal for this script)
## ------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", ask = FALSE)
}

pacman::p_load(ape, dplyr, tidyr, ggplot2)

## ------------------------------------------------------------
## Load custom helper functions
## ------------------------------------------------------------
functions_file <- "/enadisk/tempor/merlat/exogaptwo/bin/exostats_functions.R"
if (file.exists(functions_file)) {
  source(functions_file)
} else {
  stop("Custom functions file not found: ", functions_file)
}

## ------------------------------------------------------------
## Input / output paths
## ------------------------------------------------------------
newick_file <- "/enadisk/tempor/merlat/exogaptwo/cache/preprocessing/download_newick/newick/tree.nwk"
landscape_file <- "/tempor/merlat/exogaptwo/cache/repeats_annotation/merge_repeatlandscape/all/all_repeatlandscape.tsv"

out_dir <- "repeatlandscape"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## Load phylogeny and standardize tip labels
## ------------------------------------------------------------
phylo <- ape::read.tree(newick_file)
phylo$tip.label <- base::sapply(phylo$tip.label, format_genome_name)

## ------------------------------------------------------------
## Read repeat landscape table and clean numeric columns
## ------------------------------------------------------------
df <- utils::read.delim(landscape_file, header = TRUE, sep = "\t", check.names = FALSE)

df <- df %>%
  dplyr::mutate(
    div = as.integer(div),
    percentage = as.numeric(gsub(",", ".", trimws(as.character(percentage))))
  )

## Harmonize genome_id with phylogeny order (factor levels = tree tip labels)
df <- format_df(df, phylo, col = "genome_id")

## Keep only transposable elements (and optionally unknown if you want)
df <- df %>%
  dplyr::filter(type %in% c("transposable_element", "unknown"))

## Build clean_order levels
df <- factorize_by_order(df, with_TRs = FALSE)

## ------------------------------------------------------------
## Aggregate by ORDER only (genome_id x div x clean_order)
## Fill missing bins with 0 from div=0..55
## ------------------------------------------------------------
div_min <- 0L
div_max <- 55L
div_bins <- div_min:div_max

te_orders <- levels(df$clean_order)

df_order <- df %>%
  dplyr::filter(div >= div_min, div <= div_max) %>%
  dplyr::group_by(genome_id, div, clean_order) %>%
  dplyr::summarise(percentage = sum(percentage, na.rm = TRUE), .groups = "drop") %>%
  tidyr::complete(
    genome_id,
    clean_order,
    div = div_bins,
    fill = list(percentage = 0)
  ) %>%
  dplyr::mutate(
    genome_id = factor(genome_id, levels = levels(df$genome_id)),
    clean_order = factor(clean_order, levels = te_orders, exclude = NULL),
    div = as.integer(div),
    percentage = as.numeric(percentage)
  )

## ------------------------------------------------------------
## Colors per order (Unknown = darkgray)
## - RColorBrewer::Paired supports 3..12 colors
## - We assign (n_orders - 1) colors + gray for the last level
## ------------------------------------------------------------
## If you want to avoid adding a dependency, you can replace this
## block with your pre-defined colors_for_orders vector.
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer", ask = FALSE)
}
n_orders <- length(te_orders)
n_base <- max(n_orders - 1, 1)
n_paired <- min(max(n_base, 3), 12)

pal <- RColorBrewer::brewer.pal(n_paired, "Paired")
pal <- rep(pal, length.out = n_base)
pal <- c(pal, "darkgray")

colors_for_orders <- stats::setNames(pal, te_orders)
if ("Unknown" %in% te_orders) colors_for_orders["Unknown"] <- "darkgray"

## ------------------------------------------------------------
## Plot
## - stacked bars (div bins) by TE order
## - one facet per genome, ordered by phylogeny
## ------------------------------------------------------------
ncol_facets <- 4  # adjust if you want 4 per row, etc.

p <- ggplot2::ggplot(df_order, ggplot2::aes(x = div, y = percentage, fill = clean_order)) +
  ggplot2::geom_col(width = 1, color = NA) +
  ggplot2::facet_wrap(~ genome_id, ncol = ncol_facets, scales = "free_y") +
  ggplot2::scale_fill_manual(values = colors_for_orders, drop = FALSE) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
  ggplot2::scale_x_continuous(
    limits = c(div_min, div_max),
    breaks = seq(div_min, div_max - 5, by = 10),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::labs(
    x = "Divergence (%)",
    y = "Percentage",
    fill = "Order",
    title = "Repeat landscapes per genome (stacked by order)"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "italic", size = 9, hjust = 0),
    axis.text.x = ggplot2::element_text(size = 7),
    axis.text.y = ggplot2::element_text(size = 7),
    axis.title = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank()
  )

print(p)

## ------------------------------------------------------------
## Save output (single file containing all facets)
## Size is computed from number of genomes and chosen ncol
## ------------------------------------------------------------
n_genomes <- length(levels(df_order$genome_id))
nrows <- ceiling(n_genomes / ncol_facets)

width <- 10
height <- max(6, nrows * 1.2)

save_plot(
  name = "repeatlandscape",
  plot = p,
  width = width,
  height = height,
  out_dir = out_dir
)
