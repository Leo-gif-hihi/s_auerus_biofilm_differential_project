#' Microarray Analysis Script (Multi-Dataset)
#'
#' This script performs a complete analysis pipeline for microarray datasets
#' specified in input_dataset.txt. It is divided into several main phases:
#'
#' 1. Data Retrieval and Setup:
#'    - Installs and loads required Bioconductor and CRAN packages (`GEOquery`, `limma`, `affy`, `R.utils`, `ggplot2`).
#'    - Downloads raw supplementary files (TAR/CEL files) and metadata from GEO.
#'    - Unpacks the CEL files.
#'    - Generates a `targets.txt` file placeholder for manual grouping (Group and Batch).

# --- Setup Logging ---
log_file <- "Affimatrix_run.log"
cat("", file = log_file) # Initialize/Clear log file
sink(log_file, append = TRUE, split = TRUE)

# Route messages to stdout so they are also logged and split properly
message <- function(...) {
  cat(paste0(..., collapse=""), "\n", sep = "")
}
# ---------------------
#'
#' 2. Data Preprocessing and Quality Control:
#'    - Reads the microarrays using `ReadAffy` and normalizes them via `rma`.
#'    - Filters probes based on valid suffixes and `Platform_PTACC` annotations.
#'    - Filters low-expression probes by establishing a background threshold based on
#'      control probes (95th percentile), and creates a density plot (`background_distribution.png`)
#'      to visualize this distribution.
#'
#' 3. Differential Expression Analysis:
#'    - Uses `limma` to fit a linear model based on user-defined Groups and Batch.
#'    - Computes contrasts between the groups (e.g., BF5vPL5, BF7vPL7, BF9vPL9).
#'    - Annotates the results with protein accessions.
#'    - Outputs the differentially expressed genes for each contrast into CSV files
#'      (`DE_BvC.csv`, `DE_TvB.csv`, `DE_TvC.csv`).
#'
#' Note: The script expects the user to manually fill in the 'Group' and 'Batch' columns
#' in `targets.txt` after it is generated and before running the analysis phase.

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Ensure required packages are installed
required_pkgs <- c("GEOquery", "limma", "affy", "R.utils", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs)

library(GEOquery)
library(limma)
library(affy)
library(R.utils)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

# Phase 1: download data and construct targets.txt file

# Read input datasets
input_datasets <- readLines("input_dataset.txt")
input_datasets <- trimws(input_datasets)
input_datasets <- input_datasets[input_datasets != ""]

base_dir <- getwd()
annot_file_path <- file.path(base_dir, "S_aureus.na36.annot.csv")

for (dataset_name in input_datasets) {
  message(sprintf("\n\n======================================================\n=== Processing dataset: %s ===\n======================================================", dataset_name))

  dataset_dir <- file.path(base_dir, dataset_name)
  if (!dir.exists(dataset_dir)) dir.create(dataset_dir)
  setwd(dataset_dir)

  # 1. Download the supplementary files
  tar_file <- paste0(dataset_name, "_RAW.tar")
  if (!file.exists(tar_file)) {
    tryCatch(
      {
        getGEOSuppFiles(dataset_name, makeDirectory = FALSE)
      },
      error = function(e) {
        message(sprintf("Warning: Could not download supplementary files for %s: %s", dataset_name, e$message))
      }
    )
  }

  # 3. Unpack the main TAR archive
  if (file.exists(tar_file)) {
    untar(tar_file)
  }

  # 4. If files inside end in .CEL.gz, unpack them
  celfiles <- list.files(pattern = "\\.CEL\\.gz$", ignore.case = TRUE)
  if (length(celfiles) > 0) {
    sapply(celfiles, function(x) {
      out_file <- gsub("\\.gz$", "", x, ignore.case = TRUE)
      if (!file.exists(out_file)) {
        R.utils::gunzip(x, remove = FALSE)
      }
    })
  }

  # Fetch GEO dataset and metadata
  gset <- NULL
  tryCatch(
    {
      gset <- getGEO(dataset_name, AnnotGPL = TRUE, destdir = ".")
      if (length(gset) > 1) idx <- grep("GPLXXXX", attr(gset, "names")) else idx <- 1
      if (length(idx) == 0) idx <- 1
      gset <- gset[[idx]]
    },
    error = function(e) {
      message(sprintf("Warning: Could not download GEO metadata for %s: %s", dataset_name, e$message))
    }
  )

  # 5. Generate targets.txt
  if (!file.exists("targets.txt")) {
    file_names <- list.files(pattern = "\\.CEL$", ignore.case = TRUE)
    if (!is.null(gset)) {
      pd <- pData(gset)
      if ("supplementary_file_1" %in% colnames(pd)) {
        supp_files <- as.character(pd$supplementary_file_1)
        supp_file_names <- basename(supp_files)
        supp_file_names <- gsub("\\.gz$", "", supp_file_names, ignore.case = TRUE)
        # Use metadata names if they match actual files
        if (all(supp_file_names %in% file_names)) {
          file_names <- supp_file_names
        }
      }
    }

    # Just list the file names, leaving Group and Batch blank for manual entry
    targets_df <- data.frame(FileName = file_names, Group = "", Batch = "", stringsAsFactors = FALSE)
    write.table(targets_df, file = "targets.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # Phase 2: microarrays data analysis
  # 6. Read targets and process microarrays
  targets <- readTargets("targets.txt")

  if (all(targets$Group == "") || all(is.na(targets$Group))) {
    message(sprintf("Skipping analysis for %s because Group in targets.txt is empty. Please fill it and re-run.", dataset_name))
    setwd(base_dir)
    next
  }

  # Ensure contrasts.txt exists
  if (!file.exists("contrasts.txt")) {
    dummy_contrasts <- data.frame(Name = c("ConditionAvsB"), Contrast = c("ConditionA-ConditionB"), stringsAsFactors = FALSE)
    write.csv(dummy_contrasts, "contrasts.txt", row.names = FALSE, quote = FALSE)
    message(sprintf("Skipping analysis for %s: Created empty contrasts.txt. Please fill it with your contrasts and re-run.", dataset_name))
    setwd(base_dir)
    next
  }

  contrasts_df <- read.csv("contrasts.txt", stringsAsFactors = FALSE)
  if (nrow(contrasts_df) == 1 && contrasts_df$Name[1] == "ConditionAvsB") {
    message(sprintf("Skipping analysis for %s because contrasts.txt has not been customized. Please fill it and re-run.", dataset_name))
    setwd(base_dir)
    next
  }

  ab <- ReadAffy(filenames = targets$FileName)
  eset <- rma(ab)

  # Output a suffix summary of all probes
  all_probes <- rownames(eset)
  suffixes <- sub("^.*?(?=(_[a-zA-Z]_(at|st)|_(at|st))$)", "", all_probes, perl = TRUE)
  suffixes[!grepl("^(_[a-zA-Z]_(at|st)|_(at|st))$", suffixes)] <- "Other/None"
  suffix_counts <- as.data.frame(table(Suffix = suffixes))
  suffix_counts <- suffix_counts[order(-suffix_counts$Freq), ]
  write.csv(suffix_counts, paste0(dataset_name, "_probe_suffix_summary.csv"), row.names = FALSE)

  # Map probes to their extracted suffix for clean filtering
  probe_to_suffix <- setNames(suffixes, all_probes)

  # Load custom annotation from the parent directory
  f_data <- read.csv(annot_file_path, comment.char = "#", stringsAsFactors = FALSE, na.strings = "---")
  rownames(f_data) <- f_data$Probe.Set.ID

  # Format locus tags: replace /// with ;
  if ("Gene.Symbol" %in% colnames(f_data)) {
    f_data$Gene.Symbol <- as.character(f_data$Gene.Symbol)
    f_data$Gene.Symbol <- gsub("\\s*///\\s*", "; ", f_data$Gene.Symbol)
  }

  all_probes_in_eset <- rownames(eset)

  # Identify control and biological probes based on AFFX prefix
  control_probes <- all_probes_in_eset[grepl("^AFFX", all_probes_in_eset, ignore.case = TRUE)]
  biological_probes <- setdiff(all_probes_in_eset, control_probes)

  # Filter biological probes by allowed Tier 1 and Tier 2 suffixes
  allowed_suffixes <- c("_at", "_a_at")
  valid_probes <- biological_probes[probe_to_suffix[biological_probes] %in% allowed_suffixes]

  if (length(control_probes) > 0) {
    # Identify positive and negative controls based on probe IDs
    is_pos_control <- grepl("dap|lys|phe|thr", control_probes, ignore.case = TRUE)
    pos_control_probes <- control_probes[is_pos_control]
    neg_control_probes <- control_probes[!is_pos_control]

    # Use negative control probes for background if available, else all controls
    bg_ref_probes <- if (length(neg_control_probes) > 0) neg_control_probes else control_probes

    # Calculate detection p-values using MAS5 algorithm to determine expression
    message("Calculating MAS5 detection p-values for independent filtering...")
    calls <- mas5calls(ab)
    det_pvals <- assayData(calls)[["se.exprs"]]

    # Determine if valid_probes are expressed (p-value < 0.05)
    exprs_valid_pvals <- det_pvals[valid_probes, , drop = FALSE]
    is_expressed <- exprs_valid_pvals < 0.05

    # Determine minimum group size for filtering (default to 1 if unknown)
    min_samples <- 1
    if ("Group" %in% colnames(targets) && any(targets$Group != "")) {
      min_samples <- min(table(targets$Group[targets$Group != ""]))
      if (min_samples < 1) min_samples <- 1
    }

    passed_counts <- rowSums(is_expressed, na.rm = TRUE)

    # Keep probes expressed (p < 0.05) in at least `min_samples` samples
    expressed_probes <- valid_probes[passed_counts >= min_samples]
    filtered_out_probes <- valid_probes[passed_counts < min_samples]

    # Draw plot to visualize control data and split biological data
    pos_control_exprs <- as.vector(exprs(eset)[pos_control_probes, , drop = FALSE])
    neg_control_exprs <- as.vector(exprs(eset)[neg_control_probes, , drop = FALSE])
    passed_exprs <- as.vector(exprs(eset)[expressed_probes, , drop = FALSE])
    filtered_exprs <- as.vector(exprs(eset)[filtered_out_probes, , drop = FALSE])

    plot_data <- data.frame(
      Expression = c(pos_control_exprs, neg_control_exprs, passed_exprs, filtered_exprs),
      Probe_Type = factor(c(
        rep("Positive Control", length(pos_control_exprs)),
        rep("Negative Control", length(neg_control_exprs)),
        rep("Biological (Passed)", length(passed_exprs)),
        rep("Biological (Filtered Out)", length(filtered_exprs))
      ))
    )

    p <- ggplot(plot_data, aes(x = Expression, fill = Probe_Type)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(
        title = "Expression Distribution: Control vs Biological Probes",
        subtitle = paste("Passed:", length(expressed_probes), "| Filtered Out:", length(filtered_out_probes)),
        x = "Expression Level",
        y = "Density",
        fill = "Probe Type"
      ) +
      scale_fill_manual(values = c(
        "Positive Control" = "#00BFC4",
        "Negative Control" = "#F8766D",
        "Biological (Passed)" = "#7CAE00",
        "Biological (Filtered Out)" = "#C77CFF"
      ))

    ggsave(paste0(dataset_name, "_background_distribution.png"), plot = p, width = 8, height = 6, dpi = 300)

    # Save probes that failed the filter to a CSV
    if (length(filtered_out_probes) > 0) {
      filtered_df <- data.frame(
        Probe_ID = filtered_out_probes,
        Passed_Samples_Count = passed_counts[filtered_out_probes]
      )

      # Add individual array p-values
      pvals_out <- as.data.frame(exprs_valid_pvals[filtered_out_probes, , drop = FALSE])
      colnames(pvals_out) <- paste0("pval_", colnames(pvals_out))
      filtered_df <- cbind(filtered_df, pvals_out)

      # Add annotations if available
      if ("Gene.Symbol" %in% colnames(f_data)) {
        filtered_df$Gene.Symbol <- f_data[filtered_out_probes, "Gene.Symbol"]
      }

      write.csv(filtered_df, paste0(dataset_name, "_filtered_out_probes.csv"), row.names = FALSE)
      message(paste("Saved", length(filtered_out_probes), "filtered out probes to filtered_out_probes.csv"))
    }

    valid_probes <- expressed_probes
  }

  # Filter out probes that don't have a locus tag (Gene.Symbol)
  if ("Gene.Symbol" %in% colnames(f_data)) {
    gene_symbols <- f_data[valid_probes, "Gene.Symbol"]
    has_locus_tag <- !is.na(gene_symbols) & trimws(gene_symbols) != "" & trimws(gene_symbols) != "---"

    num_removed <- sum(!has_locus_tag, na.rm = TRUE)
    valid_probes <- valid_probes[has_locus_tag]
    message(sprintf("Removed %d probes without locus tags. %d probes remain.", num_removed, length(valid_probes)))
  }

  eset <- eset[valid_probes, ]

  # --- Python Integration for cluster_id ---
  # 1. Export unique locus tags to a temporary CSV
  unique_locus_tags <- unique(f_data[valid_probes, "Gene.Symbol"])
  temp_map_df <- data.frame(locus_tag = unique_locus_tags, stringsAsFactors = FALSE)
  temp_file_name <- paste0(dataset_name, "_temp_gene_mapping.csv")
  write.csv(temp_map_df, temp_file_name, row.names = FALSE)

  # 2. Call the Python scripts sequentially
  message("Calling Python scripts to retrieve cluster_id mapping...")
  python_exec <- "python"

  system2(python_exec, args = c("../../add_gene_info_new_microarray.py", "--file", temp_file_name))
  system2(python_exec, args = c("../../update_protein_ids.py", "--file", temp_file_name))
  system2(python_exec, args = c("../../add_cluster_info.py", "--file", temp_file_name))

  # 3. Read back the enriched mapping
  enriched_map <- read.csv(temp_file_name, stringsAsFactors = FALSE)
  
  if ("is_pseudogene" %in% colnames(enriched_map)) {
    is_pseudo <- as.character(enriched_map$is_pseudogene)
    is_pseudo[is.na(is_pseudo)] <- "FALSE"
    num_pseudo <- sum(toupper(is_pseudo) == "TRUE")
    if (num_pseudo > 0) {
      message(sprintf("Removing %d locus tags identified as pseudogenes.", num_pseudo))
      enriched_map <- enriched_map[toupper(is_pseudo) != "TRUE", ]
    }
  }

  # 4. Map the probes to cluster_id
  probe_locus_tags <- f_data[valid_probes, "Gene.Symbol"]
  locus_to_cluster <- setNames(enriched_map$cluster_id, enriched_map$locus_tag)
  probe_cluster_ids <- locus_to_cluster[probe_locus_tags]

  # 5. Handle missing cluster_ids (Drop as per user's choice)
  has_cluster <- !is.na(probe_cluster_ids) & trimws(probe_cluster_ids) != ""
  num_no_cluster <- sum(!has_cluster, na.rm = TRUE)

  valid_probes <- valid_probes[has_cluster]
  probe_locus_tags <- probe_locus_tags[has_cluster]
  probe_cluster_ids <- probe_cluster_ids[has_cluster]
  eset <- eset[valid_probes, ]
  message(sprintf("Dropped %d probes without cluster_id. %d probes remain.", num_no_cluster, length(valid_probes)))

  # 6. Average spot ids that share the same locus_tag
  num_before <- nrow(eset)
  eset <- avereps(eset, ID = probe_locus_tags)
  num_after <- nrow(eset)
  message(sprintf("Collapsed %d duplicate probe sets using limma::avereps based on locus_tag. %d unique locus tags remain.", num_before - num_after, num_after))

  # Ensure Group is a factor with correct levels
  Batch <- factor(targets$Batch)
  Group <- factor(targets$Group, levels = unique(targets$Group))


  # -----------------DEG analysis----------------------------------
  design <- model.matrix(~ 0 + Group + Batch)
  colnames(design) <- gsub("Group", "", colnames(design))

  fit <- lmFit(eset, design)
  cont.matrix <- makeContrasts(contrasts = contrasts_df$Contrast, levels = design)
  colnames(cont.matrix) <- contrasts_df$Name
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)

  # --- Winner Selection ---
  # Find the most significant locus_tag per cluster_id
  message("Selecting most significant locus_tag for each cluster_id based on overall F-test...")
  tt_overall <- topTable(fit2, number = Inf, sort.by = "none")
  
  unique_locus_tags <- rownames(fit2)
  fit2_cluster_ids <- locus_to_cluster[unique_locus_tags]
  
  winner_df <- data.frame(
    locus_tag = unique_locus_tags,
    cluster_id = fit2_cluster_ids,
    P.Value = tt_overall$P.Value,
    stringsAsFactors = FALSE
  )
  
  winner_df <- winner_df[order(winner_df$P.Value), ]
  winner_df <- winner_df[!duplicated(winner_df$cluster_id), ]
  
  # Subset fit2 to only keep the winning locus tags
  winning_locus_tags <- winner_df$locus_tag
  fit2 <- fit2[winning_locus_tags, ]
  
  # Seamlessly swap rownames to cluster_id for all internal matrices in fit2
  new_rownames <- winner_df$cluster_id[match(rownames(fit2), winner_df$locus_tag)]
  for (name in names(fit2)) {
    if (is.matrix(fit2[[name]]) && nrow(fit2[[name]]) == length(new_rownames)) {
      rownames(fit2[[name]]) <- new_rownames
    } else if (is.vector(fit2[[name]]) && length(fit2[[name]]) == length(new_rownames)) {
      names(fit2[[name]]) <- new_rownames
    } else if (is.data.frame(fit2[[name]]) && nrow(fit2[[name]]) == length(new_rownames)) {
      rownames(fit2[[name]]) <- new_rownames
    }
  }
  message(sprintf("Selected %d winning locus tags representing unique cluster_ids.", length(winning_locus_tags)))

  # Add gene annotations to the results
  # Since we swapped to cluster_id, we use enriched_map for annotation
  unique_enriched <- enriched_map[!duplicated(enriched_map$cluster_id) & !is.na(enriched_map$cluster_id) & trimws(enriched_map$cluster_id) != "", ]
  rownames(unique_enriched) <- unique_enriched$cluster_id

  fit2$genes <- unique_enriched[rownames(fit2), , drop = FALSE]

  # Output differentially expressed genes
  vis_dir <- "visualizations"
  if (!dir.exists(vis_dir)) dir.create(vis_dir)
  # dataset_name is provided by the loop

  # --- PCA Plots (Before and After Batch Correction) ---

  generate_pca_plot <- function(expr_matrix, title, filename) {
    cat(paste("Generating", title, "...\n"))
    pca_res <- prcomp(t(expr_matrix), scale. = TRUE)
    var_explained <- summary(pca_res)$importance[2, 1:2] * 100

    pca_df <- data.frame(
      PC1 = pca_res$x[, 1],
      PC2 = pca_res$x[, 2],
      Group = Group,
      Batch = Batch
    )

    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Batch)) +
      geom_point(size = 4, alpha = 0.8) +
      theme_minimal() +
      labs(
        title = title,
        x = sprintf("PC1 (%.1f%%)", var_explained[1]),
        y = sprintf("PC2 (%.1f%%)", var_explained[2])
      )

    ggsave(file.path(vis_dir, filename), plot = p, width = 8, height = 6, dpi = 300)
  }

  # 1. PCA Before Batch Correction
  generate_pca_plot(eset, "PCA Before Batch Correction", paste0(dataset_name, "_PCA_Before_Batch_Correction.png"))

  # 2. PCA After Batch Correction
  design_group_only <- model.matrix(~ 0 + Group)
  exprs_corrected <- removeBatchEffect(eset, batch = Batch, design = design_group_only)
  generate_pca_plot(exprs_corrected, "PCA After Batch Correction", paste0(dataset_name, "_PCA_After_Batch_Correction.png"))


  # --- Sample-to-Sample Distance Heatmaps ---
  generate_sample_distance_heatmap <- function(expr_matrix, title, filename) {
    cat(paste("Creating", title, "...\n"))
    tryCatch(
      {
        sampleDists <- dist(t(expr_matrix))
        sampleDistMatrix <- as.matrix(sampleDists)
        rownames(sampleDistMatrix) <- paste(Group, Batch, sep = "-")
        colnames(sampleDistMatrix) <- NULL

        colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

        png(file.path(vis_dir, filename), width = 800, height = 600, res = 150)
        pheatmap(sampleDistMatrix,
          clustering_distance_rows = sampleDists,
          clustering_distance_cols = sampleDists,
          col = colors,
          main = title
        )
        dev.off()
        cat("✓", title, "saved\n")
      },
      error = function(e) {
        cat("⚠ Warning: Could not create", title, ":", e$message, "\n")
      }
    )
  }

  # 1. Sample-to-Sample Distance Heatmap
  generate_sample_distance_heatmap(
    eset,
    paste("Sample-to-Sample Distances -", dataset_name),
    paste0(dataset_name, "_Sample_Distance_Heatmap.png")
  )

  # 2. Sample-to-Sample Distance Heatmap (Batch-Corrected)
  generate_sample_distance_heatmap(
    exprs_corrected,
    paste("Sample-to-Sample Distances (Batch-Corrected) -", dataset_name),
    paste0(dataset_name, "_Sample_Distance_Heatmap_Corrected.png")
  )

  # --- Expression Distribution Boxplot ---
  generate_expression_boxplot <- function(expr_matrix, title_prefix, file_prefix) {
    cat(paste("Generating", title_prefix, "expression boxplot...\n"))

    expr_vector <- as.vector(expr_matrix)
    sample_vector <- rep(colnames(expr_matrix), each = nrow(expr_matrix))
    group_vector <- rep(Group, each = nrow(expr_matrix))

    box_df <- data.frame(
      Expression = expr_vector,
      Sample = factor(sample_vector, levels = colnames(expr_matrix)),
      Group = group_vector
    )

    p_box <- ggplot(box_df, aes(x = Sample, y = Expression, fill = Group)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(title = paste(title_prefix, "- Expression Distribution"), x = "Sample", y = "Expression Level")

    ggsave(file.path(vis_dir, paste0(dataset_name, "_", file_prefix, "_Expression_Boxplot.png")), plot = p_box, width = 12, height = 6, dpi = 300)
  }

  generate_expression_boxplot(eset, "Before Batch Correction", "Before_Batch")
  generate_expression_boxplot(log2(exprs(ab)), "Raw Unnormalized (log2)", "Raw_Unnormalized")

  write_results <- function(coef_name, file_name) {
    res <- topTable(fit2, coef = coef_name, adjust = "BH", number = Inf)
    colnames(res)[colnames(res) == "logFC"] <- "log2fc"
    colnames(res)[colnames(res) == "adj.P.Val"] <- "padj"

    res_df <- as.data.frame(res)
    write.csv(res_df, file_name, row.names = TRUE)

    # A. Standard p-Value Histogram
    cat(paste("Creating p-value histogram for", coef_name, "...\n"))
    tryCatch(
      {
        png(file.path(vis_dir, paste0(dataset_name, "_PValue_Histogram_", coef_name, ".png")), width = 800, height = 600, res = 150)
        hist(res$P.Value,
          breaks = 50, col = "grey",
          main = paste("Raw p-Value Distribution -", coef_name),
          xlab = "Unadjusted p-value", ylab = "Frequency"
        )
        dev.off()
        cat("✓ p-Value histogram saved\n")
      },
      error = function(e) {
        cat("⚠ Warning: Could not create p-value histogram:", e$message, "\n")
      }
    )

    # Volcano Plot
    comparison_name <- coef_name
    res_df <- as.data.frame(res)

    # Use locus_tag if available and not empty, otherwise fallback to probe ID (rownames)
    res_df$Gene <- ifelse(!is.na(res_df$locus_tag) & trimws(res_df$locus_tag) != "",
      as.character(res_df$locus_tag),
      rownames(res_df)
    )

    # Use only the first symbol if multiple are provided (separated by ; or ,)
    res_df$Gene <- gsub("[,;].*", "", res_df$Gene)

    res_df$Significance <- "Not Significant"
    res_df$Significance[res_df$padj < 0.05 & !is.na(res_df$padj)] <- "Significant (p.adj < 0.05)"
    res_df$Significance[res_df$padj < 0.05 & abs(res_df$log2fc) > 1 & !is.na(res_df$padj)] <- "Significant & |Log2FC| > 1"

    # Create volcano plot with error handling
    tryCatch(
      {
        volcano_plot <- ggplot(res_df, aes(x = log2fc, y = -log10(padj))) +
          geom_point(aes(color = Significance), alpha = 0.6) +
          scale_color_manual(values = c(
            "Not Significant" = "grey",
            "Significant (p.adj < 0.05)" = "royalblue",
            "Significant & |Log2FC| > 1" = "red"
          )) +
          labs(
            title = paste("Volcano Plot -", dataset_name),
            subtitle = comparison_name,
            x = expression(log[2]("Fold Change")),
            y = expression(-log[10]("Adjusted p-value"))
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12)
          ) +
          geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50")

        # Add labels only if there are significant genes
        sig_genes <- subset(res_df, padj < 0.05 & !is.na(padj))
        if (nrow(sig_genes) > 0) {
          top_genes <- head(sig_genes[order(sig_genes$padj), ], 10)
          volcano_plot <- volcano_plot + geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = Inf)
        }

        # Suppress warnings when saving
        suppressWarnings(ggsave(file.path(vis_dir, paste0(dataset_name, "_Volcano_plot_", comparison_name, ".png")),
          plot = volcano_plot,
          width = 10, height = 7, dpi = 300
        ))
        cat("✓ Volcano plot saved for", comparison_name, "\n")
      },
      error = function(e) {
        cat("⚠ Warning: Could not create volcano plot for", comparison_name, ":", e$message, "\n")
      }
    )

    # MA Plot
    tryCatch(
      {
        ma_plot <- ggplot(res_df, aes(x = AveExpr, y = log2fc)) +
          geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
          scale_color_manual(values = c(
            "Not Significant" = "grey",
            "Significant (p.adj < 0.05)" = "royalblue",
            "Significant & |Log2FC| > 1" = "red"
          )) +
          labs(
            title = paste("MA Plot -", dataset_name),
            subtitle = comparison_name,
            x = "Average Expression",
            y = expression(log[2]("Fold Change"))
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12)
          ) +
          geom_hline(yintercept = 0, color = "black", linetype = "solid") +
          geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50")

        # Add labels only if there are significant genes
        if (nrow(sig_genes) > 0) {
          top_genes <- head(sig_genes[order(sig_genes$padj), ], 10)
          ma_plot <- ma_plot + geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = Inf)
        }

        suppressWarnings(ggsave(file.path(vis_dir, paste0(dataset_name, "_MA_plot_", comparison_name, ".png")),
          plot = ma_plot,
          width = 10, height = 7, dpi = 300
        ))
        cat("✓ MA plot saved for", comparison_name, "\n")
      },
      error = function(e) {
        cat("⚠ Warning: Could not create MA plot for", comparison_name, ":", e$message, "\n")
      }
    )
  }

  # --- Export Variance-Covariance Matrix ---
  cat("Extracting variance-covariance matrices for all genes...\n")

  unscaled_cov <- fit2$cov.coef
  gene_vars <- fit2$s2.post
  contrast_names <- colnames(fit2)
  num_contrasts <- length(contrast_names)

  cov_col_names <- c()
  unscaled_cov_vals <- c()

  for (i in 1:num_contrasts) {
    for (j in i:num_contrasts) {
      if (i == j) {
        cov_col_names <- c(cov_col_names, paste0("var_", contrast_names[i]))
      } else {
        cov_col_names <- c(cov_col_names, paste0("cov_", contrast_names[i], "_", contrast_names[j]))
      }
      unscaled_cov_vals <- c(unscaled_cov_vals, unscaled_cov[i, j])
    }
  }

  gene_cov_matrix <- outer(gene_vars, unscaled_cov_vals)
  colnames(gene_cov_matrix) <- cov_col_names

  # Determine cluster_id
  cluster_ids <- rownames(fit2) # Already cluster_id based on avereps

  cov_export_df <- data.frame(cluster_id = cluster_ids, stringsAsFactors = FALSE)
  cov_export_df <- cbind(cov_export_df, as.data.frame(gene_cov_matrix))

  cov_file_name <- paste0(dataset_name, "_variance_covariance.csv")
  write.csv(cov_export_df, cov_file_name, row.names = FALSE)
  cat("✓ Variance-covariance matrices saved to", cov_file_name, "\n")

  # --- Export Z-Score Covariance (Correlation) Matrix ---
  zscore_cov_matrix <- cov2cor(fit2$cov.coef)
  zscore_file_name <- paste0(dataset_name, "_zscore_covariance.csv")
  
  zscore_export_df <- as.data.frame(zscore_cov_matrix)
  write.csv(zscore_export_df, zscore_file_name, row.names = TRUE)
  cat("✓ Z-Score Covariance matrix saved to", zscore_file_name, "\n")
  # ---------------------------------------------------------

  for (c_name in contrasts_df$Name) {
    write_results(c_name, paste0(dataset_name, "_DE_", c_name, ".csv"))
  }

  # Reset working directory at the end of the loop iteration
  setwd(base_dir)
}

# Close the log file sink
sink()
