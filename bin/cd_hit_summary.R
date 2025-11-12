#!/usr/bin/env Rscript

# install.packages(c("tidyverse", "UpSetR", "VennDiagram", "patchwork", "argparse"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(UpSetR)
  library(VennDiagram)
  library(gridExtra)
  library(ggplot2)
  library(patchwork)
  library(argparse)
})

# ==========================================================
# 1. Parse the .clstr file
# ==========================================================
parse_clstr <- function(filename) {
  lines <- readLines(filename)

  data <- list()
  cluster_id <- NA

  cluster_pattern <- "^>Cluster\\s+(\\d+)"
  seq_pattern <- "^\\s*\\d+\\s+(\\d+)nt, >([^\\.]+)\\.\\.\\. (.*)$"

  for (line in lines) {
    line <- trimws(line)
    if (line == "") next

    if (grepl(cluster_pattern, line)) {
      cluster_id <- as.integer(sub(cluster_pattern, "\\1", line))
      next
    }

    if (!grepl(seq_pattern, line)) next

    length <- as.integer(sub(seq_pattern, "\\1", line))
    seq_info <- sub(seq_pattern, "\\2", line)
    end_info <- sub(seq_pattern, "\\3", line)

    # Parse source
    source_match <- str_match(seq_info, "^(RM2|HC|HLC)\\|(.+)$")
    if (!is.na(source_match[1,1])) {
      source <- source_match[1,2]
      name <- source_match[1,3]
    } else {
      source <- NA
      name <- seq_info
    }

    # Parse classification
    classification <- NA
    if (grepl("#", name)) {
      parts <- strsplit(name, "#")[[1]]
      name <- parts[1]
      classification <- parts[2]
    }

    # Representative or similarity
    if (grepl("\\*", end_info)) {
      representative <- TRUE
      similarity <- NA
    } else {
      sim_match <- str_match(end_info, "/([\\d\\.]+)%")
      similarity <- as.numeric(sim_match[1,2])
      representative <- FALSE
    }

    data[[length(data) + 1]] <- tibble(
      cluster = cluster_id,
      length = length,
      source = source,
      name = name,
      classification = classification,
      similarity = similarity,
      representative = representative
    )
  }

  bind_rows(data)
}


# ==========================================================
# 2. Add cross-tool match information
# ==========================================================
add_cross_tool_matches <- function(df) {
  cluster_sources <- df %>%
    group_by(cluster) %>%
    summarise(sources = list(unique(source)))

  df <- df %>%
    left_join(cluster_sources, by = "cluster") %>%
    mutate(has_other_tool_match = case_when(
      source == "RM2" ~ map_lgl(sources, ~ any(.x %in% c("HC", "HLC"))),
      source == "HLC" ~ map_lgl(sources, ~ "RM2" %in% .x),
      TRUE ~ NA
    )) %>%
    select(-sources)

  df
}


# ==========================================================
# 3. Filter sequences
# ==========================================================
filter_sequences <- function(df) {
  df %>%
    filter(
      (source %in% c("RM2", "HLC") & has_other_tool_match == TRUE) |
      source == "HC"
    ) %>%
    ungroup()
}


# ==========================================================
# 4. Generate all plots
# ==========================================================
generate_all_plots <- function(df, species_name, outfile) {
  # Prepare cluster sets
  rm2_classified <- df %>% filter(source == "RM2", classification != "Unknown") %>% pull(cluster) %>% unique()
  rm2_unclassified <- df %>% filter(source == "RM2", classification == "Unknown") %>% pull(cluster) %>% unique()
  hc <- df %>% filter(source == "HC") %>% pull(cluster) %>% unique()
  hlc <- df %>% filter(source == "HLC") %>% pull(cluster) %>% unique()

  # ---- UpSet plot ----
  all_clusters <- unique(c(rm2_classified, rm2_unclassified, hc, hlc))
  memberships <- lapply(all_clusters, function(cl) {
    c(if (cl %in% rm2_classified) "RM2_classified",
      if (cl %in% rm2_unclassified) "RM2_unclassified",
      if (cl %in% hc) "HiTE_confident",
      if (cl %in% hlc) "HiTE_low_confident")
  })

  upset_data <- fromList(split(all_clusters, memberships))
  upset_plot <- UpSetR::upset(
    upset_data, sets = names(upset_data),
    mainbar.y.label = "Clusters",
    order.by = "freq",
    sets.x.label = "Clusters per group"
  )

  # ---- Venn diagram ----
  rm2_clusters <- unique(df$cluster[df$source == "RM2"])
  hite_clusters <- unique(df$cluster[df$source %in% c("HC", "HLC")])
  venn.plot <- draw.pairwise.venn(
    length(rm2_clusters),
    length(hite_clusters),
    length(intersect(rm2_clusters, hite_clusters)),
    category = c("RepeatModeler2", "HiTE (HC + HLC)"),
    fill = c("skyblue", "lightgreen"),
    alpha = rep(0.7, 2)
  )

  # ---- Cluster size distribution ----
  cluster_sizes <- df %>% count(cluster)
  hist_plot <- ggplot(cluster_sizes, aes(x = n)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
    scale_y_log10() +
    labs(x = "Sequences per cluster", y = "Clusters (log scale)",
         title = "Cluster size distribution") +
    theme_minimal(base_size = 10)

  # ---- Confirmation pies ----
  df <- df %>%
    group_by(cluster) %>%
    mutate(confirmation = case_when(
      any(source == "RM2" & classification != "Unknown") ~ "confirmed by RM2 classified",
      any(source == "RM2" & classification == "Unknown") ~ "confirmed by RM2 unclassified",
      any(source == "HC") ~ "confirmed by HC",
      any(source == "HLC") ~ "confirmed by HLC",
      TRUE ~ "without confirmation"
    )) %>%
    ungroup()

  subsets <- list(
    "RepeatModeler classified" = df %>% filter(source == "RM2", classification != "Unknown"),
    "RepeatModeler unclassified" = df %>% filter(source == "RM2", classification == "Unknown"),
    "HiTE confident" = df %>% filter(source == "HC"),
    "HiTE low confident" = df %>% filter(source == "HLC")
  )

  colors <- c(
    "without confirmation" = "lightgray",
    "confirmed by RM2 classified" = "#1f77b4",
    "confirmed by RM2 unclassified" = "#aec7e8",
    "confirmed by HC" = "#2ca02c",
    "confirmed by HLC" = "#98df8a"
  )

  pie_plots <- map(subsets, function(subset) {
    counts <- subset %>%
      count(confirmation) %>%
      complete(confirmation = names(colors), fill = list(n = 0))
    ggplot(counts, aes(x = "", y = n, fill = confirmation)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      scale_fill_manual(values = colors) +
      theme_void() +
      theme(legend.position = "none")
  })

  pie_grid <- wrap_plots(pie_plots, ncol = 2) +
    plot_annotation(title = "Confirmation types by subset")

  # ---- Assemble all ----
  ggsave(outfile, pie_grid / hist_plot, width = 8.27, height = 11.69, dpi = 300)
}


# ==========================================================
# 5. Main execution
# ==========================================================
parser <- ArgumentParser(description = "CD-HIT clustering summary in R (A4 report)")
parser$add_argument("-i", "--input", required=TRUE, help="Input .clstr file")
parser$add_argument("-n", "--name", required=TRUE, help="Species name")
parser$add_argument("-o", "--output", required=TRUE, help="Output prefix")
args <- parser$parse_args()

df <- parse_clstr(args$input)
df <- add_cross_tool_matches(df)
df_filtered <- filter_sequences(df)

write_tsv(df, paste0(args$output, "_parsed.tsv"))
write_tsv(df_filtered, paste0(args$output, "_filtered.tsv"))

generate_all_plots(df, args$name, paste0(args$output, "_summary_A4.png"))

cat("✅ Analysis complete.\n")
cat(paste0(" - Parsed table: ", args$output, "_parsed.tsv\n"))
cat(paste0(" - Filtered table: ", args$output, "_filtered.tsv\n"))
cat(paste0(" - Summary A4 report: ", args$output, "_summary_A4.png\n"))
