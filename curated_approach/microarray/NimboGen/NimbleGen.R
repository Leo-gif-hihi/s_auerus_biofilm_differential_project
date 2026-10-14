#' Microarray Analysis Script (Multi-Dataset, NimbleGen)
#'
#' This script performs a complete analysis pipeline for microarray datasets
#' specified in input_dataset.txt. It is designed for NimbleGen arrays (.pair files).
#'
#' 1. Data Retrieval and Setup:
#'    - Installs and loads required Bioconductor and CRAN packages (`GEOquery`, `limma`, `R.utils`, `ggplot2`).
#'    - Downloads raw supplementary files and metadata from GEO.
#'    - Unpacks the .pair.gz files.
#'    - Generates a `targets.txt` file placeholder for manual grouping (Group and Batch).

# --- Setup Logging ---
log_file <- "NimbleGen_run.log"
cat("", file = log_file) # Initialize/Clear log file
sink(log_file, append = TRUE, split = TRUE)

# Route messages to stdout so they are also logged and split properly
message <- function(...) {
  cat(paste0(..., collapse=""), "\n", sep = "")
}
# ---------------------
#'
#' 2. Data Preprocessing and Quality Control:
#'    - Reads the microarrays using `limma::read.maimages` (generic/NimbleGen mode).
#'    - Normalizes using `backgroundCorrect` (normexp) and `normalizeBetweenArrays` (quantile).
#'    - Filters low-expression probes by establishing a background threshold based on
#'      the 95th percentile of the control probes (RANDOM).
#'
#' 3. Differential Expression Analysis:
#'    - Uses `limma` to fit a linear model based on user-defined Groups and Batch.
#'    - Computes contrasts between the groups (e.g., BF5vPL5).
#'    - Outputs the differentially expressed genes for each contrast into CSV files.
#'

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Ensure required packages are installed
required_pkgs <- c("GEOquery", "limma", "R.utils", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs)

library(GEOquery)
library(limma)
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

  # 4. If files inside end in .pair.gz, unpack them
  pairfiles <- list.files(pattern = "\\.pair\\.gz$", ignore.case = TRUE)
  if (length(pairfiles) > 0) {
    sapply(pairfiles, function(x) {
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
    file_names <- list.files(pattern = "\\.pair$", ignore.case = TRUE)
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

  message("Reading NimbleGen .pair files...")
  # Read single-channel NimbleGen pair files
  raw_data <- read.maimages(targets$FileName,
    source = "generic", columns = list(E = "PM"),
    annotation = c("PROBE_ID", "SEQ_ID", "GENE_EXPR_OPTION")
  )

  message("Normalizing data (backgroundCorrect + quantile)...")
  eset_bg <- backgroundCorrect(raw_data, method = "normexp")
  eset_probe <- normalizeBetweenArrays(eset_bg, method = "quantile")

  message("Identifying RANDOM control probes for background thresholding...")
  is_control_probe <- grepl("RANDOM", eset_probe$genes$SEQ_ID, ignore.case = TRUE) | grepl("RANDOM", eset_probe$genes$GENE_EXPR_OPTION, ignore.case = TRUE)
  control_exprs_probe <- eset_probe$E[is_control_probe, , drop = FALSE]

  # Calculate background threshold per array (95th percentile of control probes)
  bg_thresholds <- apply(control_exprs_probe, 2, function(x) quantile(x, probs = 0.95, na.rm = TRUE))

  message("Averaging replicate probes into gene-level expression...")
  eset <- avereps(eset_probe, ID = eset_probe$genes$SEQ_ID)
  all_features_in_eset <- rownames(eset$E)

  # Use annotation embedded in the NimbleGen files
  # The SEQ_ID column in NimbleGen files corresponds to the locus tag / Gene Symbol
  message("Using internal annotation from NimbleGen files...")
  f_data <- data.frame(
    Probe.Set.ID = rownames(eset$E),
    Gene.Symbol = eset$genes$SEQ_ID,
    Gene.Title = eset$genes$SEQ_ID, # Fallback, as title isn't always in pair files
    stringsAsFactors = FALSE
  )
  rownames(f_data) <- f_data$Probe.Set.ID

  # Identify control and biological features in the AVERAGED dataset
  is_control <- grepl("RANDOM", eset$genes$SEQ_ID, ignore.case = TRUE) | grepl("RANDOM", eset$genes$GENE_EXPR_OPTION, ignore.case = TRUE)

  biological_features <- all_features_in_eset[!is_control]
  valid_probes <- biological_features

  if (length(bg_thresholds) > 0 && length(biological_features) > 0) {
    message(sprintf("Found %d control probes (yielding %d summarized control features) and %d summarized biological genes.", sum(is_control_probe), sum(is_control), length(biological_features)))

    # Check if each biological feature is above the threshold in each sample
    bio_exprs <- eset$E[biological_features, , drop = FALSE]
    is_expressed <- sweep(bio_exprs, 2, bg_thresholds, ">")

    # Determine minimum group size for filtering (default to 1 if unknown)
    min_samples <- 1
    if ("Group" %in% colnames(targets) && any(targets$Group != "")) {
      min_samples <- min(table(targets$Group[targets$Group != ""]))
      if (min_samples < 1) min_samples <- 1
    }

    passed_counts <- rowSums(is_expressed, na.rm = TRUE)

    # Keep features expressed in at least `min_samples` samples
    expressed_features <- biological_features[passed_counts >= min_samples]
    filtered_out_features <- biological_features[passed_counts < min_samples]

    # Draw plot to visualize control data and split biological data
    passed_exprs <- as.vector(bio_exprs[expressed_features, , drop = FALSE])
    filtered_exprs <- as.vector(bio_exprs[filtered_out_features, , drop = FALSE])
    ctrl_exprs_vec <- as.vector(control_exprs_probe)

    plot_data <- data.frame(
      Expression = c(ctrl_exprs_vec, passed_exprs, filtered_exprs),
      Probe_Type = factor(c(
        rep("Control (RANDOM)", length(ctrl_exprs_vec)),
        rep("Biological (Passed)", length(passed_exprs)),
        rep("Biological (Filtered Out)", length(filtered_exprs))
      ))
    )

    p <- ggplot(plot_data, aes(x = Expression, fill = Probe_Type)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(
        title = "Expression Distribution: Control vs Biological Probes",
        subtitle = paste("Passed:", length(expressed_features), "| Filtered Out:", length(filtered_out_features)),
        x = "Expression Level (Normalized)",
        y = "Density",
        fill = "Probe Type"
      ) +
      scale_fill_manual(values = c(
        "Control (RANDOM)" = "#F8766D",
        "Biological (Passed)" = "#7CAE00",
        "Biological (Filtered Out)" = "#C77CFF"
      ))

    ggsave(paste0(dataset_name, "_background_distribution.png"), plot = p, width = 8, height = 6, dpi = 300)

    # Save features that failed the filter to a CSV
    if (length(filtered_out_features) > 0) {
      filtered_df <- data.frame(
        Locus_Tag = filtered_out_features,
        Passed_Samples_Count = passed_counts[filtered_out_features]
      )

      status_out <- as.data.frame(is_expressed[filtered_out_features, , drop = FALSE])
      colnames(status_out) <- paste0("passed_in_", colnames(status_out))
      filtered_df <- cbind(filtered_df, status_out)

      # Add annotations if available
      if ("Gene.Symbol" %in% colnames(f_data)) {
        gs <- as.character(f_data[filtered_out_features, "Gene.Symbol"])
        gs <- gsub("\\s*///\\s*", "; ", gs)
        filtered_df$Gene.Symbol <- gs
      }

      write.csv(filtered_df, paste0(dataset_name, "_filtered_out_features.csv"), row.names = FALSE)
      message(paste("Saved", length(filtered_out_features), "filtered out features to filtered_out_features.csv"))
    }

    valid_probes <- expressed_features
  }
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
  generate_pca_plot(eset$E, "PCA Before Batch Correction", paste0(dataset_name, "_PCA_Before_Batch_Correction.png"))

  # 2. PCA After Batch Correction
  design_group_only <- model.matrix(~ 0 + Group)
  exprs_corrected <- removeBatchEffect(eset$E, batch = Batch, design = design_group_only)
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
    eset$E,
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

  generate_expression_boxplot(eset$E, "Before Batch Correction", "Before_Batch")
  generate_expression_boxplot(log2(raw_data$E), "Raw Unnormalized (log2)", "Raw_Unnormalized")

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
