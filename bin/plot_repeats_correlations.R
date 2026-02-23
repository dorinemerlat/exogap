#!/usr/bin/env Rscript

## ------------------------------------------------------------
## Correlation plots per subtree of a phylogeny
## ------------------------------------------------------------

# To delete
setwd("/tempor/merlat/exogaptwo")

## ------------------------------------------------------------
## Packages
## ------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", ask = FALSE)
}
pacman::p_load(ape, phylogram, dplyr, ggplot2, ggtext, glue)

## ------------------------------------------------------------
## Load custom functions
## ------------------------------------------------------------
functions_file <- "/enadisk/tempor/merlat/exogaptwo/bin/exostats_functions.R"
if (file.exists(functions_file)) {
  source(functions_file)
} else {
  stop("Custom functions file not found: ", functions_file)
}

## ------------------------------------------------------------
## Input files
## ------------------------------------------------------------
newick_file <- "/enadisk/tempor/merlat/exogaptwo/cache/preprocessing/download_newick/newick/tree.nwk"
stats_file  <- "/tempor/merlat/exogaptwo/cache/preprocessing/merge_stats/all/all.stats"
summary_all_file <- "/enadisk/tempor/merlat/exogaptwo/cache/repeats_annotation/merge_summaries/repeats_summary_all_only_repetitive_elements/repeats_summary_all_only_repetitive_elements.tsv"

## ------------------------------------------------------------
## Load phylogeny and standardize genome names
## ------------------------------------------------------------
phylo <- ape::read.tree(newick_file)
phylo$tip.label <- sapply(phylo$tip.label, format_genome_name)

## ------------------------------------------------------------
## Load data tables and harmonize with phylogeny tip labels
## ------------------------------------------------------------
stats <- read_summary(stats_file, phylo)

summary_all <- read_summary(summary_all_file, phylo) %>%
  dplyr::left_join(stats, by = "genome_id") %>%
  dplyr::rename(size = sum_len)

## Scale numeric columns for readable axis units
res <- select_unit(summary_all, list(size = "pb", load = "n", N50 = "pb"))
summary_all <- res$df
units <- res$units

## ------------------------------------------------------------
## Function: draw one correlation plot per subtree
## ------------------------------------------------------------
draw_correlation <- function(df, x_column, y_column, tree, units, out_dir) {
  
  ## ape::subtrees() returns all clades/subtrees in the phylogeny
  subtrees_list <- ape::subtrees(ape::as.phylo(tree))
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (subtree in subtrees_list) {
    
    ## Parent node label (clade name) + children (tips)
    parent <- subtree[["node.label"]][[1]]
    childs <- subtree[["tip.label"]]
    
    ## Clean labels (remove extra suffix after '|', if any)
    parent <- sub("\\|.*$", "", parent)
    childs <- sub("\\|.*$", "", childs)
    
    ## Skip tiny subtrees (not enough points to correlate)
    if (length(childs) <= 3) next
    
    ## Subset dataframe to genomes present in this subtree
    sub_df <- df[df$genome_id %in% childs, , drop = FALSE]
    
    ## Skip if requested columns are missing (robust for future changes)
    if (!x_column %in% colnames(sub_df)) next
    if (!y_column %in% colnames(sub_df)) next
    
    ## Compute correlations (suppress warnings for ties / small samples)
    spearman <- suppressWarnings(stats::cor.test(
      sub_df[[x_column]], sub_df[[y_column]],
      method = "spearman", exact = FALSE
    ))
    
    pearson <- suppressWarnings(stats::cor.test(
      sub_df[[x_column]], sub_df[[y_column]],
      method = "pearson"
    ))
    
    ## Add correlation text label to the dataframe
    sub_df <- dplyr::mutate(
      sub_df,
      label = glue::glue(
        "*r* = {round(pearson$estimate, 3)}, *p* = {round(pearson$p.value, 4)}\n\n",
        "*&rho;* = {round(spearman$estimate, 3)}, *p* = {round(spearman$p.value, 4)}"
      )
    )
    
    ## Axis labels
    x_label <- ifelse(x_column == "coverage", "Coverage (%)", glue::glue("Load{units$load}"))
    y_label <- ifelse(y_column == "size", glue::glue("Assembly size{units$size}"), "N50")
    
    ## Build plot
    p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data[[x_column]], y = .data[[y_column]])) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
      ggplot2::ggtitle(parent) +
      ggplot2::ylab(y_label) +
      ggplot2::xlab(x_label) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(size = 0.5, color = "grey80"),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(size = 1, color = "black"),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      ) +
      ggtext::geom_richtext(
        data = sub_df,
        ggplot2::aes(label = label, x = -Inf, y = Inf),
        hjust = -0.1,
        vjust = 1.1,
        fill = scales::alpha("white", 0.5)
      )
    
    ## Output filename
    plot_name <- glue::glue("correlation_{parent}_{tolower(x_column)}_vs_{tolower(y_column)}")
    
    ## Save outputs (pdf + svg via your custom save_plot)
    save_plot(plot_name, p, width = 4, height = 4, out_dir = out_dir)
  }
}

## ------------------------------------------------------------
## Run correlations for multiple y-columns
## ------------------------------------------------------------
run_correlations <- function(df, phylo, units, y_columns = c("size", "N50")) {
  for (y in y_columns) {
    out_dir <- glue::glue("correlation_{tolower(y)}")
    draw_correlation(df, "coverage", y, phylo, units, out_dir)
    draw_correlation(df, "load",     y, phylo, units, out_dir)
  }
}

## ------------------------------------------------------------
## Main
## ------------------------------------------------------------
run_correlations(df = summary_all, phylo = phylo, units = units)

