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
log_file <- "Agilent_run.log"
cat("", file = log_file) # Initialize/Clear log file
sink(log_file, append = TRUE, split = TRUE)

# Route messages to stdout so they are also logged and split properly
message <- function(...) {
  cat(paste0(..., collapse=""), "\n", sep = "")
}
# ---------------------
#'
#' 2. Data Preprocessing and Quality Control:
#'    - Reads the microarrays using `limma::read.maimages`.
#'    - Applies `normexp` background correction and `loess` within-array normalization.
#'    - Filters out control probes and unexpressed biological probes using an intensity threshold.
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
#' Note: The script expects the user to manually fill in the 'Cy3' column
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
library(R.utils)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

# Phase 1: download data and construct targets.txt file

# Read input datasets
input_datasets <- readLines("input_dataset.txt", warn = FALSE)
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

  # 4. If files inside end in .txt.gz, unpack them
  txtfiles <- list.files(pattern = "\\.txt\\.gz$", ignore.case = TRUE)
  if (length(txtfiles) > 0) {
    sapply(txtfiles, function(x) {
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
    file_names <- list.files(pattern = "\\.txt$", ignore.case = TRUE)
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

    # Generate targets with Cy5 = "gDNA" and Cy3 left blank for manual input
    targets_df <- data.frame(FileName = file_names, Cy3 = "", Cy5 = "gDNA", Batch = "", stringsAsFactors = FALSE)
    write.table(targets_df, file = "targets.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # Phase 2: microarrays data analysis
  # 6. Read targets and process microarrays
  targets <- readTargets("targets.txt")

  if (all(targets$Cy3 == "") || all(is.na(targets$Cy3))) {
    message(sprintf("Skipping analysis for %s because Cy3 in targets.txt is empty. Please fill it and re-run.", dataset_name))
    setwd(base_dir)
    next
  }

  # Check if files have mismatched line counts
  line_counts <- sapply(targets$FileName, function(f) {
    if (file.exists(f)) R.utils::countLines(f)[1] else 0
  })

  if (length(unique(line_counts)) > 1) {
    message("Warning: Files have different row counts. Performing inner join on ProbeName...")

    common_probes <- NULL
    for (f in targets$FileName) {
      if (!file.exists(f)) next
      lines <- readLines(f, n = 100)
      header_idx <- grep("^FEATURES", lines)[1]
      d <- read.delim(f, skip = header_idx - 1, stringsAsFactors = FALSE)
      if (is.null(common_probes)) {
        common_probes <- d$ProbeName
      } else {
        common_probes <- intersect(common_probes, d$ProbeName)
      }
    }

    temp_dir <- file.path(tempdir(), paste0("truncated_arrays_", dataset_name))
    dir.create(temp_dir, showWarnings = FALSE)
    new_filenames <- file.path(temp_dir, basename(targets$FileName))

    for (i in seq_along(targets$FileName)) {
      f <- targets$FileName[i]
      if (!file.exists(f)) next
      lines <- readLines(f)
      header_idx <- grep("^FEATURES", lines)[1]
      header <- lines[1:header_idx]

      d <- read.delim(f, skip = header_idx - 1, stringsAsFactors = FALSE)

      d_filtered <- d[d$ProbeName %in% common_probes, ]
      # Sort to ensure identical order
      d_filtered <- d_filtered[match(common_probes, d_filtered$ProbeName), ]

      out_f <- new_filenames[i]
      writeLines(header, out_f)
      write.table(d_filtered, out_f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    targets$FileName <- new_filenames
    message(sprintf("Found %d common probes. Created joined arrays in: %s", length(common_probes), temp_dir))
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

  RG <- read.maimages(targets, source = "agilent", green.only = FALSE)
  RG <- backgroundCorrect(RG, method = "normexp", offset = 50)
  MA <- normalizeWithinArrays(RG, method = "loess")

  # Filter control probes and unexpressed biological probes
  all_probes <- MA$genes$ProbeName

  # ControlType in Agilent: 1 = positive, -1 = negative, 0 = biological
  if ("ControlType" %in% colnames(MA$genes)) {
    is_pos_control <- !is.na(MA$genes$ControlType) & MA$genes$ControlType == 1
    is_neg_control <- !is.na(MA$genes$ControlType) & MA$genes$ControlType == -1
    is_biological <- !is.na(MA$genes$ControlType) & MA$genes$ControlType == 0
  } else {
    is_pos_control <- rep(FALSE, nrow(MA))
    is_neg_control <- rep(FALSE, nrow(MA))
    is_biological <- rep(TRUE, nrow(MA))
  }

  pos_control_probes <- all_probes[is_pos_control]
  neg_control_probes <- all_probes[is_neg_control]
  biological_probes <- all_probes[is_biological]

  # Use negative control probes for background threshold using A values (average log2 intensity)
  A_values <- MA$A

  # We find the 95th percentile of the negative control probes A values as background
  neg_control_A <- as.vector(A_values[is_neg_control, , drop = FALSE])
  if (length(neg_control_A) > 0 && sum(!is.na(neg_control_A)) > 0) {
    bg_threshold <- quantile(neg_control_A, 0.95, na.rm = TRUE)
  } else {
    bg_threshold <- 5 # Fallback
  }

  message(sprintf("Calculated background expression threshold (A value): %.2f", bg_threshold))

  # Determine if biological probes are expressed
  is_expressed <- A_values > bg_threshold
  is_expressed[is.na(is_expressed)] <- FALSE
  exprs_valid_flags <- is_expressed[is_biological, , drop = FALSE]

  min_samples <- 1
  if ("Cy3" %in% colnames(targets) && any(targets$Cy3 != "")) {
    min_samples <- min(table(targets$Cy3[targets$Cy3 != ""]))
    if (min_samples < 1) min_samples <- 1
  }

  passed_counts <- rowSums(exprs_valid_flags, na.rm = TRUE)
  expressed_biological <- biological_probes[passed_counts >= min_samples]
  filtered_out_biological <- biological_probes[passed_counts < min_samples]

  # Draw plot to visualize control data and split biological data
  pos_control_exprs <- as.vector(A_values[is_pos_control, , drop = FALSE])
  passed_exprs <- as.vector(A_values[match(expressed_biological, all_probes), , drop = FALSE])
  filtered_exprs <- as.vector(A_values[match(filtered_out_biological, all_probes), , drop = FALSE])

  plot_data <- data.frame(
    Expression = c(pos_control_exprs, neg_control_A, passed_exprs, filtered_exprs),
    Probe_Type = factor(c(
      rep("Positive Control", length(pos_control_exprs)),
      rep("Negative Control", length(neg_control_A)),
      rep("Biological (Passed)", length(passed_exprs)),
      rep("Biological (Filtered Out)", length(filtered_exprs))
    ))
  )

  p <- ggplot(plot_data, aes(x = Expression, fill = Probe_Type)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(
      title = "Expression Distribution (A Values): Control vs Biological Probes",
      subtitle = paste("Passed:", length(expressed_biological), "| Filtered Out:", length(filtered_out_biological)),
      x = "Average Expression (A value)",
      y = "Density",
      fill = "Probe Type"
    ) +
    scale_fill_manual(values = c(
      "Positive Control" = "#00BFC4",
      "Negative Control" = "#F8766D",
      "Biological (Passed)" = "#7CAE00",
      "Biological (Filtered Out)" = "#C77CFF"
    )) +
    geom_vline(xintercept = bg_threshold, linetype = "dashed", color = "red")

  ggsave(paste0(dataset_name, "_background_distribution.png"), plot = p, width = 8, height = 6, dpi = 300)

  if (length(filtered_out_biological) > 0) {
    filtered_df <- data.frame(
      ProbeName = filtered_out_biological,
      Passed_Samples_Count = passed_counts[filtered_out_biological]
    )

    if ("SystematicName" %in% colnames(MA$genes)) {
      sys_names <- MA$genes$SystematicName[match(filtered_out_biological, MA$genes$ProbeName)]
      filtered_df$SystematicName <- sys_names
    }

    write.csv(filtered_df, paste0(dataset_name, "_filtered_out_probes.csv"), row.names = FALSE)
    message(paste("Saved", length(filtered_out_biological), "filtered out probes to filtered_out_probes.csv"))
  }

  # Keep only expressed biological probes
  valid_indices <- match(expressed_biological, all_probes)

  # Combine SystematicName and GeneName and extract only standard locus tags
  has_sys <- "SystematicName" %in% colnames(MA$genes)
  has_gn <- "GeneName" %in% colnames(MA$genes)

  if (has_sys || has_gn) {
    clean_locus_tags <- function(sys_names) {
      sapply(sys_names, function(s) {
        if (is.na(s) || trimws(s) == "" || s == "---") {
          return(NA)
        }
        # Remove quotes that might prevent regex matching
        s <- gsub("[\"']", "", s)

        # Split by both semicolon and comma
        parts <- strsplit(s, "[;,]")[[1]]
        parts <- trimws(parts)
        parts <- parts[parts != "" & parts != "*"]

        valid_parts <- parts[grepl("^[A-Za-z0-9-]+_[A-Za-z0-9]+$", parts)]
        valid_parts <- valid_parts[!grepl("^(IO|TilingRegion)", valid_parts, ignore.case = TRUE)]

        valid_parts <- unique(valid_parts)

        if (length(valid_parts) > 0) {
          return(paste(valid_parts, collapse = ", "))
        } else {
          return(NA)
        }
      }, USE.NAMES = FALSE)
    }

    sys_col <- if (has_sys) MA$genes$SystematicName else rep(NA, nrow(MA))
    gn_col <- if (has_gn) MA$genes$GeneName else rep(NA, nrow(MA))

    # Combine both columns to ensure we catch all locus tags, then store back into SystematicName
    combined_names <- paste(sys_col, gn_col, sep = ";")
    MA$genes$SystematicName <- clean_locus_tags(combined_names)
  }

  # Filter out probes that don't have a locus tag (SystematicName)
  if ("SystematicName" %in% colnames(MA$genes)) {
    gene_symbols <- MA$genes$SystematicName[valid_indices]
    has_locus_tag <- !is.na(gene_symbols) & trimws(gene_symbols) != "" & trimws(gene_symbols) != "---"

    num_removed <- sum(!has_locus_tag, na.rm = TRUE)
    valid_indices <- valid_indices[has_locus_tag]
    message(sprintf("Removed %d probes without locus tags. %d probes remain.", num_removed, length(valid_indices)))
  }

  MA <- MA[valid_indices, ]

  # --- Python Integration for cluster_id ---
  # 1. Export unique locus tags to a temporary CSV
  # Use SystematicName if available, otherwise it falls back
  gene_symbols <- if ("SystematicName" %in% colnames(MA$genes)) MA$genes$SystematicName[valid_indices] else rep(NA, length(valid_indices))
  unique_locus_tags <- unique(gene_symbols[!is.na(gene_symbols)])
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
  locus_to_cluster <- setNames(enriched_map$cluster_id, enriched_map$locus_tag)
  probe_cluster_ids <- locus_to_cluster[gene_symbols]

  # 5. Handle missing cluster_ids (Drop as per user's choice)
  has_cluster <- !is.na(probe_cluster_ids) & trimws(probe_cluster_ids) != ""
  num_no_cluster <- sum(!has_cluster, na.rm = TRUE)

  valid_indices <- valid_indices[has_cluster]
  gene_symbols <- gene_symbols[has_cluster]
  probe_cluster_ids <- probe_cluster_ids[has_cluster]
  MA <- MA[valid_indices, ]
  message(sprintf("Dropped %d probes without cluster_id. %d probes remain.", num_no_cluster, length(valid_indices)))

  # 6. Average spot ids that share the same locus_tag
  num_before <- nrow(MA)
  MA <- avereps(MA, ID = gene_symbols)
  num_after <- nrow(MA)
  message(sprintf("Collapsed %d duplicate probes using limma::avereps based on locus_tag. %d unique locus tags remain.", num_before - num_after, num_after))

  # Ensure Group is a factor with correct levels
  Batch <- factor(targets$Batch)

  # Sanitize Cy3 and Cy5 names to valid R variable names (e.g. starting with numbers is invalid)
  original_cy3 <- as.character(targets$Cy3)
  targets$Cy3 <- make.names(original_cy3)
  if ("Cy5" %in% colnames(targets)) {
    targets$Cy5 <- make.names(as.character(targets$Cy5))
  }
  Group <- factor(targets$Cy3, levels = unique(targets$Cy3))

  # Update user contrast strings to match the sanitized names
  unique_orig <- unique(original_cy3[original_cy3 != ""])
  # Sort by length descending to prevent substring replacement bugs
  unique_orig <- unique_orig[order(nchar(unique_orig), decreasing = TRUE)]

  for (orig in unique_orig) {
    valid <- make.names(orig)
    if (orig != valid) {
      contrasts_df$Contrast <- gsub(orig, valid, contrasts_df$Contrast, fixed = TRUE)
    }
  }

  # -----------------DEG analysis----------------------------------
  # For two color array with common reference, create design matrix
  design <- modelMatrix(targets, ref = "gDNA")

  fit <- lmFit(MA, design)
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
  generate_pca_plot(MA$M, "PCA Before Batch Correction", paste0(dataset_name, "_PCA_Before_Batch_Correction.png"))

  # 2. PCA After Batch Correction
  design_group_only <- model.matrix(~ 0 + Group)
  exprs_corrected <- removeBatchEffect(MA$M, batch = Batch, design = design_group_only)
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
    MA$M,
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

  generate_expression_boxplot(MA$M, "Normalized M-values", "Normalized_M_Values")
  generate_expression_boxplot(log2(RG$R / RG$G), "Raw Unnormalized (log2 R/G)", "Raw_Unnormalized")

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

    # Use locus_tag if available and not empty, otherwise fallback to cluster_id (rownames)
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
