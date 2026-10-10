#!/usr/bin/env Rscript

# generate_custom_forest_plot.R
# This script reads the DESeq2 contrast results and dataset coefficients
# and generates a forest plot for a custom list of target genes.

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-g", "--genes"),
    type = "character", default = "purH,icaA,pbp2,icaD,putP,gltB",
    help = "Comma-separated list of target genes (Representative_Gene or cluster_id) e.g., 'agrA,sarA,YP_498609.1'", metavar = "character"
  ),
  make_option(c("-f", "--file"),
    type = "character", default = NULL,
    help = "Text file with list of target genes (one per line)", metavar = "character"
  ),
  make_option(c("-r", "--results"),
    type = "character", default = "results_sample_info-biofilm-to-biofilm-deseq2only/DESeq2_contrast_results.csv",
    help = "Path to DESeq2_contrast_results.csv [default= %default]", metavar = "character"
  ),
  make_option(c("-c", "--coefs"),
    type = "character", default = "results_sample_info-biofilm-to-biofilm-deseq2only/dataset_coefficients.csv",
    help = "Path to dataset_coefficients.csv [default= %default]", metavar = "character"
  ),
  make_option(c("-m", "--model"),
    type = "character", default = "FE",
    help = "Meta-analysis model to plot: 3LRE, RE, FE, GAMMAGEE, or DIRIGNORELRE [default= %default]", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "custom_forest_plot.png",
    help = "Output filename for the plot [default= %default]", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$genes) && is.null(opt$file)) {
  print_help(opt_parser)
  stop("You must provide target genes using either --genes or --file", call. = FALSE)
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# Load data
if (!file.exists(opt$results) || !file.exists(opt$coefs)) {
  stop("Could not find results or coefs files. Please check the --results and --coefs paths.\nDefault looks for 'DESeq2_contrast_results.csv' in the current directory.", call. = FALSE)
}

cat("Loading data...\n")
res_df <- read.csv(opt$results, stringsAsFactors = FALSE)
coef_df <- read.csv(opt$coefs, stringsAsFactors = FALSE)

# Get target genes
target_genes <- c()
if (!is.null(opt$genes)) {
  target_genes <- unlist(strsplit(opt$genes, ","))
}
if (!is.null(opt$file)) {
  if (file.exists(opt$file)) {
    lines <- readLines(opt$file)
    target_genes <- unique(c(target_genes, lines[lines != ""]))
  } else {
    stop(paste("Gene list file not found:", opt$file), call. = FALSE)
  }
}
target_genes <- trimws(target_genes)
target_genes <- target_genes[target_genes != ""]

# Find matching genes (match against cluster_id or Representative_Gene)
matched_res <- res_df %>%
  filter(cluster_id %in% target_genes | Representative_Gene %in% target_genes)

if (nrow(matched_res) == 0) {
  stop("None of the target genes were found in the results file. Check spelling or ID types.", call. = FALSE)
}

cat("Found", nrow(matched_res), "matching genes out of", length(target_genes), "requested.\n")

model_type <- toupper(opt$model)
if (!model_type %in% c("3LRE", "RE", "FE", "GAMMAGEE", "DIRIGNORELRE")) {
  stop("Model must be one of: 3LRE, RE, FE, GAMMAGEE, DIRIGNORELRE", call. = FALSE)
}

lfc_col_meta <- paste0("Meta_", model_type, "_LFC")
se_col_meta <- paste0("Meta_", model_type, "_SE")
overall_label <- paste0("Overall (", model_type, ")")

if (!lfc_col_meta %in% colnames(res_df)) {
  stop(paste("Column", lfc_col_meta, "not found in results file."), call. = FALSE)
}

# Identify comparison columns in coef_df
coef_cols <- colnames(coef_df)
lfc_cols <- coef_cols[grepl("_LFC_adj$", coef_cols)]
comp_names <- sub("_LFC_adj$", "", lfc_cols)

plot_data_list <- list()

for (i in seq_len(nrow(matched_res))) {
  gene_id <- matched_res$cluster_id[i]
  gene_name <- matched_res$Representative_Gene[i]
  gene_label <- paste0(gene_name, " (", gene_id, ")")

  # Individual comparisons
  gene_coefs <- coef_df %>% filter(cluster_id == gene_id)

  for (comp in comp_names) {
    comp_lfc_col <- paste0(comp, "_LFC_adj")
    comp_se_col <- paste0(comp, "_SE")

    if (comp_lfc_col %in% colnames(gene_coefs) && comp_se_col %in% colnames(gene_coefs)) {
      lfc_val <- gene_coefs[[comp_lfc_col]][1]
      se_val <- gene_coefs[[comp_se_col]][1]

      if (!is.na(lfc_val) && !is.na(se_val)) {
        plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
          Gene = gene_label,
          Comparison = comp,
          Type = "Individual",
          LFC = lfc_val,
          SE = se_val,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Meta-analysis result
  meta_lfc <- matched_res[[lfc_col_meta]][i]
  meta_se <- matched_res[[se_col_meta]][i]

  if (!is.na(meta_lfc) && !is.na(meta_se)) {
    plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
      Gene = gene_label,
      Comparison = overall_label,
      Type = "Meta-Analysis",
      LFC = meta_lfc,
      SE = meta_se,
      stringsAsFactors = FALSE
    )
  }
}

plot_data <- do.call(rbind, plot_data_list)

if (is.null(plot_data) || nrow(plot_data) == 0) {
  stop("No valid LFC/SE data found to plot for the requested genes.", call. = FALSE)
}

# Fix ordering (Overall at the bottom)
comps_unique <- setdiff(unique(plot_data$Comparison), overall_label)
plot_data$Comparison <- factor(plot_data$Comparison, levels = c(rev(comps_unique), overall_label))

ordered_genes <- unique(plot_data$Gene)
plot_data$Gene <- factor(plot_data$Gene, levels = ordered_genes)

# Plot
cat("Generating plot...\n")
p_forest <- ggplot(plot_data, aes(x = LFC, y = Comparison, color = Type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(shape = Type, size = Type)) +
  geom_errorbarh(aes(xmin = LFC - 1.96 * SE, xmax = LFC + 1.96 * SE), height = 0.2) +
  scale_color_manual(values = c("Individual" = "black", "Meta-Analysis" = "red")) +
  scale_shape_manual(values = c("Individual" = 16, "Meta-Analysis" = 18)) +
  scale_size_manual(values = c("Individual" = 2, "Meta-Analysis" = 3)) +
  facet_wrap(~Gene, scales = "free_y", ncol = min(3, length(ordered_genes))) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = paste("Custom Forest Plot of Log2 Fold Changes (", model_type, ")", sep = ""),
    x = "Log2 Fold Change (\u00B1 95% CI)", y = ""
  )

# Calculate height based on number of genes/comparisons
n_genes <- length(ordered_genes)
n_rows <- ceiling(n_genes / 3)
# Base height is around 4 per row, adjust based on number of comparisons if there are many
n_comps <- length(comps_unique)
height_per_row <- max(4, n_comps * 0.25)
plot_height <- max(6, n_rows * height_per_row)

ggsave(opt$output, p_forest, width = 14, height = plot_height, dpi = 300)
cat("✓ Forest Plot successfully saved to", opt$output, "\n")
