# DESeq2 Pipeline - Independent Models Meta-Analysis Approach
#
# This script performs differential expression analysis using DESeq2
# independently on all GSE datasets. It then calculates a meta-analysis
# effect size across datasets using metafor.
#
# INPUT FILES:
# - annot.csv (metadata indicating biofilm direction)
# - raw/*.xlsx or raw/*.csv (raw count data)
# - sample_info/*.xlsx or sample_info/*.csv (sample metadata for each dataset)

# --- 0. Parse Arguments & Setup Output Directory ---
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org")
}
library(optparse)

option_list <- list(
  make_option(c("-s", "--sample_dir"),
    type = "character", default = "sample_info-biofilm-to-planktonic",
    help = "Directory containing sample info files [default= %default]", metavar = "character"
  ),
  make_option(c("-m", "--min_datasets"),
    type = "integer", default = 0,
    help = "Minimum number of datasets a gene must appear in to be analyzed. 0 means all datasets. [default= %default]", metavar = "integer"
  ),
  make_option(c("--use_cwb"),
    action = "store_true", default = FALSE,
    help = "Enable Cluster Wild Bootstrapping for the 3LRE model (computationally heavy) [default= %default]"
  ),
  make_option(c("--cwb_r"),
    type = "integer", default = 399,
    help = "Number of bootstrap replications if --use_cwb is enabled [default= %default]", metavar = "integer"
  ),
  make_option(c("--cwb_pval_thresh"),
    type = "numeric", default = 0.05,
    help = "Only run CWB on genes with an initial 3LRE p-value below this threshold [default= %default]", metavar = "numeric"
  ),
  make_option(c("--use_cwb_stouffer"),
    action = "store_true", default = FALSE,
    help = "Enable Cluster Wild Bootstrapping for the RE Stouffer model (computationally heavy) [default= %default]"
  ),
  make_option(c("--cwb_stouffer_r"),
    type = "integer", default = 399,
    help = "Number of bootstrap replications if --use_cwb_stouffer is enabled [default= %default]", metavar = "integer"
  ),
  make_option(c("--cwb_stouffer_pval_thresh"),
    type = "numeric", default = 0.05,
    help = "Only run CWB on genes with an initial RE Stouffer p-value below this threshold [default= %default]", metavar = "numeric"
  ),
  make_option(c("--use_pb"),
    action = "store_true", default = FALSE,
    help = "Enable Parametric Bootstrapping for the Gamma GLMM model (computationally heavy) [default= %default]"
  ),
  make_option(c("--pb_r"),
    type = "integer", default = 500,
    help = "Number of bootstrap replications if --use_pb is enabled [default= %default]", metavar = "integer"
  ),
  make_option(c("--pb_pval_thresh"),
    type = "numeric", default = 0.05,
    help = "Only run PB on genes with an initial GLMM p-value below this threshold [default= %default]", metavar = "numeric"
  ),
  make_option(c("--ora_padj"),
    type = "numeric", default = 0.05,
    help = "Adjusted p-value threshold for ORA gene selection [default= %default]", metavar = "numeric"
  ),
  make_option(c("--ora_lfc"),
    type = "numeric", default = 0.5,
    help = "Absolute log2FC threshold for ORA gene selection [default= %default]", metavar = "numeric"
  ),
  make_option(c("--gaf_file"),
    type = "character", default = "QuickGO-annotations-1764241851574-20251127.gaf",
    help = "Path to GAF file for GO annotations [default= %default]", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

sample_dir <- opt$sample_dir
sample_dir_name <- basename(sample_dir)

out_dir <- paste0("results_", sample_dir_name)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# --- 0.5 Initialize Logging ---
sink(file.path(out_dir, "DESeq2_independent_analysis_log.txt"), split = TRUE)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("DESEQ2 ANALYSIS - INDEPENDENT MODELS APPROACH\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# --- 1. Installation and Loading ---
cat("Installing/loading required packages...\n")

install_if_needed <- function(package, bioc = FALSE) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(package)
    } else {
      install.packages(package, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_needed("DESeq2", bioc = TRUE)
install_if_needed("limma", bioc = TRUE)
install_if_needed("pheatmap")
install_if_needed("RColorBrewer")
install_if_needed("openxlsx")
install_if_needed("ggplot2")
install_if_needed("ggrepel")
install_if_needed("dplyr")
install_if_needed("tidyr")
install_if_needed("stringr")
install_if_needed("apeglm", bioc = TRUE)
install_if_needed("ashr")
install_if_needed("metafor")
install_if_needed("clubSandwich")
install_if_needed("EmpiricalBrownsMethod", bioc = TRUE)
install_if_needed("glmmTMB")

library(DESeq2)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(apeglm)
library(ashr)
library(metafor)
library(clubSandwich)
library(EmpiricalBrownsMethod)
library(glmmTMB)

# --- Create Visualization Directory ---
vis_dir <- file.path(out_dir, "Visualizations_Independent")
if (!dir.exists(vis_dir)) {
  dir.create(vis_dir, recursive = TRUE)
  cat("Created visualization directory:", vis_dir, "\n")
}

# --- 2. Load Annotation Rules ---
cat("\nLoading annotation rules from annot.csv...\n")
annot <- read.csv("annot.csv", stringsAsFactors = FALSE)
annot_map <- annot %>%
  select(GSE.ID, Modified.in.biofilm.formation.compare.to.control..description.) %>%
  rename(GSE = GSE.ID, Direction = 2)
cat("Annotation rules loaded successfully.\n")

# --- 3. Load and Combine Data ---
cat("\n=== LOADING AND MERGING DATASETS ===\n")
raw_files <- list.files("raw", full.names = TRUE, pattern = "\\.(csv|xlsx)$")
sample_info_files <- list.files(sample_dir, full.names = TRUE, pattern = "\\.(csv|xlsx)$")
if (length(sample_info_files) == 0) {
  stop(paste("CRITICAL ERROR: No .csv or .xlsx files found in", sample_dir))
}

# Check for duplicate sample info files with the same GSE ID
gse_ids_all <- str_extract(basename(sample_info_files), "GSE[0-9]+")
dup_gse <- unique(gse_ids_all[!is.na(gse_ids_all) & duplicated(gse_ids_all)])
if (length(dup_gse) > 0) {
  dup_files <- sample_info_files[gse_ids_all %in% dup_gse]
  stop(paste0(
    "CRITICAL ERROR: Multiple sample info files found for the same GSE ID(s): ",
    paste(dup_gse, collapse = ", "), "\nConflicting files:\n  ",
    paste(basename(dup_files), collapse = "\n  ")
  ))
}

all_counts <- list()
all_coldata <- list()
all_gene_mappings <- list()

for (si_file in sample_info_files) {
  gse_id <- str_extract(basename(si_file), "GSE[0-9]+")
  cat("\nProcessing", gse_id, "...\n")

  if (grepl("\\.xlsx$", si_file)) {
    si <- read.xlsx(si_file)
  } else {
    si <- read.csv(si_file, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # Find matching raw file
  raw_file <- raw_files[grepl(gse_id, raw_files)]
  if (length(raw_file) == 0) {
    cat("  ⚠ No raw data found for", gse_id, "- skipping\n")
    next
  }
  if (length(raw_file) > 1) {
    stop(paste0(
      "CRITICAL ERROR: Multiple raw data files found for ", gse_id, ":\n  ",
      paste(basename(raw_file), collapse = "\n  ")
    ))
  }

  if (grepl("\\.xlsx$", raw_file)) {
    raw <- read.xlsx(raw_file)
  } else {
    raw <- read.csv(raw_file, stringsAsFactors = FALSE, check.names = FALSE)
  }

  if (!"cluster_id" %in% colnames(raw)) {
    stop(paste0("CRITICAL ERROR: 'cluster_id' column missing in ", gse_id, " raw data (", basename(raw_file), ")"))
  }


  # Extract gene mapping for this dataset
  gene_col <- intersect(c("gene", "gene_name", "locus_tag", "gene_id", "cluster_id"), colnames(raw))[1]
  if (!is.na(gene_col)) {
    mapping_df <- raw[, c("cluster_id", gene_col)]
    colnames(mapping_df) <- c("cluster_id", "gene")
    mapping_df <- mapping_df[!is.na(mapping_df$cluster_id) & mapping_df$cluster_id != "", ]

    # Filter out rRNA locus tags (e.g., SAOUHSC_R0009 or NWMN_rRNA)
    is_rrna <- grepl("^(gene-)?SAOUHSC_R[0-9]+|NWMN_rRNA", mapping_df$gene)
    if (any(is_rrna)) {
      cat("  - Removed", sum(is_rrna), "rRNA locus tags\n")
      mapping_df <- mapping_df[!is_rrna, ]
      # Filter raw to exclude these rRNA genes from counts
      raw <- raw[raw$cluster_id %in% mapping_df$cluster_id, ]
    }

    mapping_df$GSE <- gse_id
    all_gene_mappings[[gse_id]] <- mapping_df
  }

  # Identify samples present in both
  samples <- intersect(si$Sample, colnames(raw))
  if (length(samples) == 0) {
    stop(paste("CRITICAL ERROR: No matching samples between info and raw data for", gse_id))
  }

  # Check for mismatches
  missing_in_raw <- setdiff(si$Sample, colnames(raw))
  missing_in_si <- setdiff(colnames(raw), si$Sample)

  if (length(missing_in_raw) > 0 || length(missing_in_si) > 0) {
    cat("  ⚠ Warning: Mismatch in samples for", gse_id, "\n")
    if (length(missing_in_raw) > 0) {
      cat("    Samples in info but missing in raw data:", paste(missing_in_raw, collapse = ", "), "\n")
    }
    if (length(missing_in_si) > 0) {
      cat("    Samples in raw data but missing in info:", paste(missing_in_si, collapse = ", "), "\n")
    }
  }

  si <- si[si$Sample %in% samples, ]

  # Map Groups and Identify Baseline/Control
  si$Original_Group <- si$Group
  si$Group <- make.names(si$Group)

  if (!gse_id %in% c("GSE162709", "GSE204884")) {
    is_control <- grepl("Control|planktonic", si$Original_Group, ignore.case = TRUE)
    if (any(is_control)) {
      si$Group[is_control] <- "Control"
    } else {
      stop(paste0("CRITICAL ERROR: No 'Control' or 'planktonic' group found for dataset ", gse_id, "."))
    }
  }

  # Batch
  if ("Batch" %in% colnames(si)) {
    si$Batch <- as.factor(si$Batch)
  } else {
    si$Batch <- as.factor("1")
  }

  si$GSE <- gse_id

  si_clean <- si[, c("Sample", "GSE", "Batch", "Group", "Original_Group")]

  # Make sample names globally unique
  unique_samples <- paste0(gse_id, "_", samples)
  si_clean$Sample <- unique_samples

  # Extract counts
  counts <- raw[, c("cluster_id", samples)]
  colnames(counts) <- c("cluster_id", unique_samples)
  counts[, unique_samples] <- apply(counts[, unique_samples, drop = FALSE], 2, function(x) as.numeric(as.character(x)))

  # Aggregate by cluster_id
  counts_agg <- counts %>%
    group_by(cluster_id) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(n_datasets = 1)

  all_counts[[gse_id]] <- counts_agg
  all_coldata[[gse_id]] <- si_clean

  cat("  ✓ Processed", gse_id, "successfully. (", length(samples), "samples, ", nrow(counts_agg), "unique genes)\n")
}

cat("\n=== PROCESSING MICROARRAY DATASETS FOR GENE COUNT ===\n")
microarray_list_file <- file.path(sample_dir, "microarray.txt")
microarray_dir_map <- list()
microarray_nt_map <- list()
microarray_nc_map <- list()
microarray_shared_map <- list()
microarray_files <- character()

if (file.exists(microarray_list_file)) {
  ma_lines <- readLines(microarray_list_file)
  ma_lines <- ma_lines[ma_lines != ""]

  if (length(ma_lines) > 0) {
    allowed_microarray <- character()
    for (l in ma_lines) {
      parts <- strsplit(l, "\t")[[1]]
      file_name <- trimws(parts[1])
      allowed_microarray <- c(allowed_microarray, file_name)
      if (length(parts) >= 2) microarray_dir_map[[file_name]] <- trimws(parts[2])
      if (length(parts) >= 3) microarray_nt_map[[file_name]] <- as.numeric(trimws(parts[3]))
      if (length(parts) >= 4) microarray_nc_map[[file_name]] <- as.numeric(trimws(parts[4]))
      if (length(parts) >= 5) microarray_shared_map[[file_name]] <- trimws(parts[5])
    }

    cand_files <- file.path("microarray", allowed_microarray)
    microarray_files <- cand_files[file.exists(cand_files)]
  }
}
microarray_gses <- unique(str_extract(basename(microarray_files), "GSE[0-9]+"))
microarray_gses <- microarray_gses[!is.na(microarray_gses)]
if (length(microarray_gses) > 0) {
  for (gse_id in microarray_gses) {
    gse_files <- microarray_files[grepl(gse_id, microarray_files)]

    gse_clusters <- c()
    for (f in gse_files) {
      df <- tryCatch(
        {
          read.csv(f, stringsAsFactors = FALSE)
        },
        error = function(e) NULL
      )
      if (!is.null(df) && "cluster_id" %in% colnames(df)) {
        gene_col <- intersect(c("gene", "gene_name", "locus_tag", "gene_id"), colnames(df))[1]
        if (!is.na(gene_col)) {
          mapping_df <- df[, c("cluster_id", gene_col)]
          colnames(mapping_df) <- c("cluster_id", "gene")
          mapping_df <- mapping_df[!is.na(mapping_df$cluster_id) & mapping_df$cluster_id != "", ]

          is_rrna <- grepl("^(gene-)?SAOUHSC_R[0-9]+|NWMN_rRNA", mapping_df$gene)
          if (any(is_rrna)) {
            mapping_df <- mapping_df[!is_rrna, ]
            df <- df[!df$cluster_id %in% df$cluster_id[is_rrna], ]
          }

          mapping_df$GSE <- gse_id
          all_gene_mappings[[paste0(gse_id, "_", basename(f))]] <- mapping_df
        }

        gse_clusters <- unique(c(gse_clusters, df$cluster_id))
      }
    }

    gse_clusters <- gse_clusters[!is.na(gse_clusters) & gse_clusters != ""]

    if (length(gse_clusters) > 0) {
      counts_agg <- data.frame(cluster_id = gse_clusters, n_datasets = 1, stringsAsFactors = FALSE)
      all_counts[[gse_id]] <- counts_agg
      cat("  ✓ Processed microarray", gse_id, "successfully. (", length(gse_clusters), "unique genes)\n")
    } else {
      cat("  ⚠ No valid cluster_ids found for microarray", gse_id, "\n")
    }
  }
}

cat("\n=== PROCESSING GENE MAPPINGS ===\n")
gene_map_combined <- bind_rows(all_gene_mappings)
gene_map_combined <- gene_map_combined %>%
  filter(!is.na(gene) & gene != "") %>%
  mutate(
    gene_clean = gsub("^gene-", "", gene),
    is_locus = nchar(gene_clean) > 6 & grepl("_", gene_clean)
  )

cluster_annotation <- gene_map_combined %>%
  group_by(cluster_id) %>%
  summarise(
    All_Genes_Locus_Tags = paste(unique(gene_clean), collapse = "/"),
    Representative_Gene = {
      symbols <- unique(gene_clean[!is_locus])
      locus_tags <- unique(gene_clean[is_locus])
      if (length(symbols) > 0) {
        paste(symbols, collapse = "/")
      } else if (length(locus_tags) > 0) {
        locus_tags[1]
      } else {
        cluster_id[1]
      }
    },
    .groups = "drop"
  )

cat("\n=== MERGING ALL DATASETS ===\n")
combined_counts_df <- bind_rows(all_counts)
combined_counts_df <- combined_counts_df %>%
  group_by(cluster_id) %>%
  summarise(
    n_datasets = sum(n_datasets, na.rm = TRUE),
    across(everything(), ~ sum(.x, na.rm = TRUE))
  )

# Filter for Core/Minimum Genome
expected_datasets <- length(all_counts)
min_ds <- opt$min_datasets
if (min_ds == 0) min_ds <- expected_datasets

cat("\n=== GENE FILTERING ===\n")
cat("Filtering for genes present in at least", min_ds, "datasets.\n")

valid_genes <- combined_counts_df$n_datasets >= min_ds
cat("Genes passing filter:", sum(valid_genes), "\n")

if (sum(valid_genes) == 0) {
  stop("CRITICAL ERROR: No genes passed the dataset occurrence filter!")
}

combined_counts_df <- combined_counts_df[valid_genes, ]
combined_counts_df$n_datasets <- NULL
combined_counts_df <- as.data.frame(combined_counts_df)

valid_rows <- !is.na(combined_counts_df$cluster_id) & combined_counts_df$cluster_id != ""
combined_counts_df <- combined_counts_df[valid_rows, , drop = FALSE]

rownames(combined_counts_df) <- make.names(combined_counts_df$cluster_id, unique = TRUE)
combined_counts_df$cluster_id <- NULL

cts <- as.matrix(combined_counts_df)
cts <- apply(cts, 2, function(x) round(as.numeric(x)))
rownames(cts) <- rownames(combined_counts_df)

combined_coldata <- bind_rows(all_coldata)
rownames(combined_coldata) <- combined_coldata$Sample

common_samples <- intersect(colnames(cts), rownames(combined_coldata))
cts <- cts[, common_samples]
combined_coldata <- combined_coldata[common_samples, ]

# Prepare colData factors
colData <- combined_coldata
colData$GSE <- as.factor(colData$GSE)

sign_map <- list()

# --- 4. DESeq2 Analysis (Independent per Dataset) ---
cat("\n=== RUNNING DESEQ2 ANALYSIS INDEPENDENTLY ===\n")

unique_gses <- unique(as.character(colData$GSE))
n_genes <- nrow(cts)

# We will collect LFC, SE, and p-values as columns dynamically
lfc_list <- list()
se_list <- list()
pval_list <- list()
neff_list <- list()
nt_list <- list()
nc_list <- list()
comparison_gse_map <- c() # To map comparison name to GSE ID/cluster for meta-analysis

dataset_dds_list <- list()

for (gse_id in unique_gses) {
  cat("\nRunning DESeq2 for dataset:", gse_id, "\n")

  # Subset colData and cts
  gse_colData <- colData[colData$GSE == gse_id, ]

  gse_colData$Group <- as.factor(gse_colData$Group)
  if ("Control" %in% levels(gse_colData$Group)) {
    gse_colData$Group <- relevel(gse_colData$Group, ref = "Control")
  } else if (!gse_id %in% c("GSE162709", "GSE204884")) {
    cat("  ⚠ 'Control' not found in Group for", gse_id, "- skipping\n")
    next
  }

  gse_colData$Batch <- as.factor(gse_colData$Batch)
  gse_colData$Group <- droplevels(gse_colData$Group)
  gse_colData$Batch <- droplevels(gse_colData$Batch)

  gse_counts <- cts[, rownames(gse_colData), drop = FALSE]

  # Remove genes with 0 counts in this specific dataset
  keep_genes <- rowSums(gse_counts) > 0
  gse_counts_filtered <- gse_counts[keep_genes, , drop = FALSE]

  cat("  Genes with > 0 counts:", nrow(gse_counts_filtered), "\n")

  if (nrow(gse_counts_filtered) == 0) {
    cat("  ⚠ No genes to analyze. Skipping.\n")
    next
  }

  # Determine formula
  if (length(unique(gse_colData$Batch)) > 1) {
    design_formula <- ~ Batch + Group
    cat("  Design:", as.character(design_formula), "\n")
  } else {
    design_formula <- ~Group
    cat("  Design:", as.character(design_formula), "\n")
  }

  # Create DDS
  dds_gse <- DESeqDataSetFromMatrix(
    countData = gse_counts_filtered,
    colData = gse_colData,
    design = design_formula
  )

  # Filter low counts for this dataset
  min_group_size <- min(table(gse_colData$Group)[table(gse_colData$Group) > 0])
  if (is.infinite(min_group_size)) min_group_size <- 1
  keep <- rowSums(counts(dds_gse) >= 10) >= min_group_size
  dds_gse <- dds_gse[keep, ]
  cat("  Genes after low-count filtering:", nrow(dds_gse), "\n")

  # Run DESeq2
  dds_gse <- tryCatch(
    {
      DESeq(dds_gse, quiet = TRUE)
    },
    error = function(e) {
      cat("  ⚠ Error in DESeq2 for", gse_id, ":", e$message, "\n")
      return(NULL)
    }
  )

  if (is.null(dds_gse)) next
  dataset_dds_list[[gse_id]] <- dds_gse

  # Determine custom comparisons based on dataset
  if (gse_id == "GSE204884") {
    # The groups in gse_colData are syntactically valid names.
    comparisons <- list(
      list(comp = "cody.3h", ref = "Control.3h"),
      list(comp = "ccpa.3h", ref = "Control.3h"),
      list(comp = "codyccpa.3h", ref = "Control.3h"),
      list(comp = "cody.6h", ref = "Control.6h"),
      list(comp = "ccpa.6h", ref = "Control.6h"),
      list(comp = "codyccpa.6h", ref = "Control.6h")
    )
    cat("  Using custom comparisons for GSE204884 (time-matched controls)\n")
  } else if (gse_id == "GSE162709") {
    comparisons <- list(
      list(comp = "Biofilm.SH1000.TSB", ref = "Control.SH1000.TSB"),
      list(comp = "Biofilm.SH1000.BHI", ref = "Control.SH1000.BHI"),
      list(comp = "Biofilm.USA300.TSB", ref = "Control.USA300.TSB"),
      list(comp = "Biofilm.USA300.BHI", ref = "Control.USA300.BHI")
    )
    cat("  Using custom comparisons for GSE162709 (strain and media-matched controls)\n")
  } else {
    # Default: compare all other groups to "Control"
    ref_group <- "Control"
    comp_groups <- setdiff(levels(dds_gse$Group), ref_group)
    comparisons <- lapply(comp_groups, function(cg) list(comp = cg, ref = ref_group))
  }

  dir_info <- annot_map$Direction[annot_map$GSE == gse_id]

  for (comp_info in comparisons) {
    comp_group <- comp_info$comp
    ref_group <- comp_info$ref

    # Check if groups exist in this dataset
    if (!comp_group %in% levels(dds_gse$Group) || !ref_group %in% levels(dds_gse$Group)) {
      cat("  ⚠ Groups", comp_group, "or", ref_group, "not found. Skipping.\n")
      next
    }

    comparison_name <- paste0(gse_id, "_", comp_group, "_vs_", ref_group)
    cat("  Extracting comparison:", comparison_name, "\n")

    res_gse <- tryCatch(
      {
        results(dds_gse, contrast = c("Group", comp_group, ref_group))
      },
      error = function(e) {
        cat("  ⚠ Error extracting contrast:", e$message, "\n")
        return(NULL)
      }
    )

    if (is.null(res_gse)) next

    # Determine sign
    sign <- 0
    if (gse_id == "GSE204884") {
      if (grepl("cody(?!ccpa)", comp_group, ignore.case = TRUE, perl = TRUE)) {
        sign <- 1
      } else {
        sign <- -1
      }
    } else if (length(dir_info) > 0 && !is.na(dir_info) && dir_info != "") {
      if (grepl("increased", dir_info, ignore.case = TRUE)) {
        sign <- 1
      } else if (grepl("decreased", dir_info, ignore.case = TRUE)) sign <- -1
    }

    if (sign == 0) {
      cat("  ⚠ Could not determine sign for", comparison_name, ". Setting sign to NA.\n")
      sign <- NA
    } else {
      cat("  Applying sign multiplier:", sign, "\n")
    }
    sign_map[[comparison_name]] <- sign

    # Store LFC, SE, and P-value
    genes_in_res <- rownames(res_gse)

    lfc_vec <- rep(NA, n_genes)
    names(lfc_vec) <- rownames(cts)
    lfc_vec[genes_in_res] <- res_gse$log2FoldChange * sign

    se_vec <- rep(NA, n_genes)
    names(se_vec) <- rownames(cts)
    se_vec[genes_in_res] <- res_gse$lfcSE

    pval_vec <- rep(NA, n_genes)
    names(pval_vec) <- rownames(cts)
    pval_vec[genes_in_res] <- res_gse$pvalue

    # Calculate N_eff
    n_t <- sum(gse_colData$Group == comp_group)
    n_c <- sum(gse_colData$Group == ref_group)
    n_eff <- (n_t * n_c) / (n_t + n_c)

    lfc_list[[comparison_name]] <- lfc_vec
    se_list[[comparison_name]] <- se_vec
    pval_list[[comparison_name]] <- pval_vec
    neff_list[[comparison_name]] <- n_eff
    nt_list[[comparison_name]] <- n_t
    nc_list[[comparison_name]] <- n_c
    comparison_gse_map[comparison_name] <- gse_id
  }

  # --- Extract Variance-Covariance Matrix for DESeq2 ---
  cat("  Extracting variance-covariance matrices for", gse_id, "...\n")
  build_contrast <- function(factor_name, num_level, den_level, res_names) {
    factor_terms <- grep(paste0("^", factor_name, "_"), res_names, value = TRUE)
    if (length(factor_terms) > 0) {
      base_level <- sub(".*_vs_", "", factor_terms[1])
      if (num_level == base_level) {
        c_num <- rep(0, length(res_names))
      } else {
        num_name <- paste0(factor_name, "_", num_level, "_vs_", base_level)
        c_num <- as.numeric(res_names == num_name)
      }
      if (den_level == base_level) {
        c_den <- rep(0, length(res_names))
      } else {
        den_name <- paste0(factor_name, "_", den_level, "_vs_", base_level)
        c_den <- as.numeric(res_names == den_name)
      }
      return(c_num - c_den)
    } else {
      # Fallback if no specific factor terms found
      return(rep(0, length(res_names)))
    }
  }

  res_names <- resultsNames(dds_gse)
  C_mat <- matrix(0, nrow = length(res_names), ncol = length(comparisons))
  comp_names_gse <- sapply(comparisons, function(x) paste0(gse_id, "_", x$comp, "_vs_", x$ref))
  colnames(C_mat) <- comp_names_gse

  for (k in seq_along(comparisons)) {
    C_mat[, k] <- build_contrast("Group", comparisons[[k]]$comp, comparisons[[k]]$ref, res_names)
  }
  out_vcov_dir <- file.path(out_dir, "deseq2_var_cov_results")
  vcov_file <- file.path(out_vcov_dir, paste0(gse_id, "_variance_covariance.csv"))

  cov_mats_gse <- list()
  valid_genes <- rownames(dds_gse)

  if (file.exists(vcov_file)) {
    cat("  Found existing variance-covariance matrix file for", gse_id, ". Loading it...
")
    df_vcov <- read.csv(vcov_file, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(df_vcov))) {
      gene <- as.character(df_vcov$cluster_id[i])
      mat <- matrix(NA, nrow = length(comp_names_gse), ncol = length(comp_names_gse), dimnames = list(comp_names_gse, comp_names_gse))
      for (k in seq_along(comp_names_gse)) {
        comp_name <- comp_names_gse[k]
        var_col <- paste0("var_", comp_name)
        if (var_col %in% colnames(df_vcov)) mat[k, k] <- df_vcov[i, var_col]
        if (k < length(comp_names_gse)) {
          for (l in (k + 1):length(comp_names_gse)) {
            comp_name2 <- comp_names_gse[l]
            cov_col <- paste0("cov_", comp_name, "_", comp_name2)
            if (cov_col %in% colnames(df_vcov)) {
              mat[k, l] <- df_vcov[i, cov_col]
              mat[l, k] <- df_vcov[i, cov_col]
            }
          }
        }
      }
      cov_mats_gse[[gene]] <- mat
    }
  } else {
    for (i in seq_along(valid_genes)) {
      gene_cov <- tryCatch(
        {
          V <- DESeq2:::covarianceMatrix(dds_gse, rowNumber = i)
          t(C_mat) %*% V %*% C_mat
        },
        error = function(e) {
          matrix(NA, nrow = ncol(C_mat), ncol = ncol(C_mat), dimnames = list(comp_names_gse, comp_names_gse))
        }
      )
      cov_mats_gse[[valid_genes[i]]] <- gene_cov
    }

    df_vcov <- data.frame(cluster_id = valid_genes, stringsAsFactors = FALSE)
    for (k in seq_along(comp_names_gse)) {
      comp_name <- comp_names_gse[k]
      df_vcov[[paste0("var_", comp_name)]] <- sapply(cov_mats_gse, function(x) x[k, k])
      if (k < length(comp_names_gse)) {
        for (l in (k + 1):length(comp_names_gse)) {
          comp_name2 <- comp_names_gse[l]
          df_vcov[[paste0("cov_", comp_name, "_", comp_name2)]] <- sapply(cov_mats_gse, function(x) x[k, l])
        }
      }
    }

    if (!dir.exists(out_vcov_dir)) dir.create(out_vcov_dir, recursive = TRUE)
    write.csv(df_vcov, vcov_file, row.names = FALSE)
  }

  # --- Extract Z-Score Covariance (Correlation) Matrix for DESeq2 ---
  zcov_mats_gse <- list()
  for (i in seq_along(valid_genes)) {
    gene_cov <- cov_mats_gse[[valid_genes[i]]]
    zcov <- tryCatch(cov2cor(gene_cov), error = function(e) {
      matrix(NA, nrow = ncol(C_mat), ncol = ncol(C_mat), dimnames = list(comp_names_gse, comp_names_gse))
    })
    zcov_mats_gse[[valid_genes[i]]] <- zcov
  }

  df_zcov <- data.frame(cluster_id = valid_genes, stringsAsFactors = FALSE)
  for (k in seq_along(comp_names_gse)) {
    comp_name <- comp_names_gse[k]
    df_zcov[[paste0("var_", comp_name)]] <- 1
    if (k < length(comp_names_gse)) {
      for (l in (k + 1):length(comp_names_gse)) {
        comp_name2 <- comp_names_gse[l]
        df_zcov[[paste0("cov_", comp_name, "_", comp_name2)]] <- sapply(zcov_mats_gse, function(x) x[k, l])
      }
    }
  }

  out_zcov_dir <- file.path(out_dir, "deseq2_zscore_cov_results")
  if (!dir.exists(out_zcov_dir)) dir.create(out_zcov_dir, recursive = TRUE)
  write.csv(df_zcov, file.path(out_zcov_dir, paste0(gse_id, "_zscore_covariance.csv")), row.names = FALSE)

  if (!exists("all_cov_mats")) all_cov_mats <- list()
  all_cov_mats[[gse_id]] <- cov_mats_gse

  if (!exists("all_zcov_mats")) all_zcov_mats <- list()
  all_zcov_mats[[gse_id]] <- zcov_mats_gse
}

cat("\n=== PROCESSING MICROARRAY DE RESULTS ===\n")
if (length(microarray_files) > 0) {
  # Pre-load microarray variance-covariance matrices
  microarray_vcov_dir <- "microarray_var_cov_results"
  ma_vcov_files <- list.files(microarray_vcov_dir, pattern = "\\.csv$", full.names = TRUE)
  ma_vcov_data <- list()
  for (f in ma_vcov_files) {
    gse_id_match <- str_extract(basename(f), "GSE[0-9]+")
    if (!is.na(gse_id_match)) {
      ma_vcov_data[[gse_id_match]] <- read.csv(f, stringsAsFactors = FALSE)
    }
  }

  # Pre-load microarray Z-score covariance matrices
  microarray_zcov_dir <- "microarry_zscore_cov_results"
  ma_zcov_files <- list.files(microarray_zcov_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!exists("all_zcov_mats")) all_zcov_mats <- list()
  for (f in ma_zcov_files) {
    gse_id_match <- str_extract(basename(f), "GSE[0-9]+")
    if (!is.na(gse_id_match)) {
      zcov_mat <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
      full_names <- paste0(gse_id_match, "_DE_", rownames(zcov_mat))
      rownames(zcov_mat) <- full_names
      colnames(zcov_mat) <- full_names
      all_zcov_mats[[gse_id_match]] <- list(dataset_level = as.matrix(zcov_mat))
    }
  }

  # Build all_cov_mats for microarray datasets based on loaded var_cov data
  cat("  Constructing variance-covariance matrices for microarray datasets...\n")
  if (!exists("all_cov_mats")) all_cov_mats <- list()
  for (v_gse_id in names(ma_vcov_data)) {
    vcov_df <- ma_vcov_data[[v_gse_id]]
    comp_cols <- grep("^var_", colnames(vcov_df), value = TRUE)
    comps <- sub("^var_", "", comp_cols) # These are suffixes, e.g. BF5vPL5
    full_comp_names <- paste0(v_gse_id, "_DE_", comps)

    cov_mats_gse <- list()
    for (i in 1:nrow(vcov_df)) {
      gene <- as.character(vcov_df$cluster_id[i])
      mat <- matrix(NA, nrow = length(comps), ncol = length(comps), dimnames = list(full_comp_names, full_comp_names))
      for (r in 1:length(comps)) {
        mat[r, r] <- vcov_df[[paste0("var_", comps[r])]][i]
        if (r < length(comps)) {
          for (c in (r + 1):length(comps)) {
            cov_val <- vcov_df[[paste0("cov_", comps[r], "_", comps[c])]][i]
            if (is.null(cov_val)) cov_val <- vcov_df[[paste0("cov_", comps[c], "_", comps[r])]][i]
            mat[r, c] <- cov_val
            mat[c, r] <- cov_val
          }
        }
      }
      cov_mats_gse[[gene]] <- mat
    }
    all_cov_mats[[v_gse_id]] <- cov_mats_gse
  }

  for (mf in microarray_files) {
    comp_name <- sub("\\.csv$", "", basename(mf)) # e.g. GSE138075_DE_BF5vPL5
    gse_id <- str_extract(comp_name, "GSE[0-9]+")

    cat("\nExtracting microarray comparison:", comp_name, "\n")
    df <- tryCatch(
      {
        read.csv(mf, stringsAsFactors = FALSE)
      },
      error = function(e) NULL
    )

    if (is.null(df) || !"cluster_id" %in% colnames(df) || !"log2fc" %in% colnames(df) || !"t" %in% colnames(df)) {
      cat("  ⚠ Missing required columns (cluster_id, log2fc, t) in", comp_name, ". Skipping.\n")
      next
    }

    # Extract SE from the exact variance matrix (diagonal)
    if (!is.null(all_cov_mats[[gse_id]])) {
      df$SE <- sapply(df$cluster_id, function(gene) {
        if (!is.null(all_cov_mats[[gse_id]][[gene]])) {
          var_val <- all_cov_mats[[gse_id]][[gene]][comp_name, comp_name]
          if (!is.na(var_val) && var_val > 0) {
            return(sqrt(var_val))
          }
        }
        return(NA)
      })
    } else {
      # Fallback to reverse-engineering SE if variance matrix is missing
      cat("  ⚠ No exact variance matrix found for", comp_name, "- falling back to log2fc/t\n")
      df$SE <- abs(df$log2fc / df$t)
      df$SE[is.infinite(df$SE) | is.nan(df$SE)] <- NA
    }

    # Determine directionality from microarray.txt
    dir_info <- microarray_dir_map[[basename(mf)]]
    sign <- 0
    if (!is.null(dir_info) && dir_info != "") {
      if (grepl("increased", dir_info, ignore.case = TRUE)) {
        sign <- 1
      } else if (grepl("decreased", dir_info, ignore.case = TRUE)) {
        sign <- -1
      }
    }

    if (sign == 0) {
      cat("  ⚠ Could not determine sign for", comp_name, ". Setting sign to NA.\n")
      sign <- NA
    } else {
      cat("  Applying sign multiplier:", sign, "\n")
    }
    sign_map[[comp_name]] <- sign

    # Extract p-value
    if (!"P.Value" %in% colnames(df)) {
      cat("  ⚠ Missing P.Value in", comp_name, "- Setting to NA.\n")
      pval_col <- rep(NA, nrow(df))
    } else {
      pval_col <- df$P.Value
    }
    df$P.Value <- pval_col

    # Calculate N_eff and shared ID
    n_t <- microarray_nt_map[[basename(mf)]]
    n_c <- microarray_nc_map[[basename(mf)]]
    if (is.null(n_t) || is.na(n_t) || is.null(n_c) || is.na(n_c)) {
      cat("  ⚠ Missing N_T or N_C for", comp_name, "- Using default N_eff = 1.5\n")
      n_eff <- 1.5
    } else {
      n_eff <- (n_t * n_c) / (n_t + n_c)
    }

    shared_id <- microarray_shared_map[[basename(mf)]]
    if (is.null(shared_id) || is.na(shared_id) || shared_id == "") {
      shared_id <- gse_id
    }

    # Enforce uniqueness instead of aggregating
    df_filtered <- df %>% filter(!is.na(cluster_id) & cluster_id %in% rownames(cts))
    if (any(duplicated(df_filtered$cluster_id))) {
      stop(paste("CRITICAL ERROR: Duplicate cluster_id found in microarray dataset", comp_name))
    }

    lfc_vec <- rep(NA, n_genes)
    names(lfc_vec) <- rownames(cts)
    lfc_vec[df_filtered$cluster_id] <- df_filtered$log2fc * sign

    se_vec <- rep(NA, n_genes)
    names(se_vec) <- rownames(cts)
    se_vec[df_filtered$cluster_id] <- df_filtered$SE

    pval_vec <- rep(NA, n_genes)
    names(pval_vec) <- rownames(cts)
    pval_vec[df_filtered$cluster_id] <- df_filtered$P.Value

    lfc_list[[comp_name]] <- lfc_vec
    se_list[[comp_name]] <- se_vec
    pval_list[[comp_name]] <- pval_vec
    neff_list[[comp_name]] <- n_eff
    nt_list[[comp_name]] <- n_t
    nc_list[[comp_name]] <- n_c
    comparison_gse_map[comp_name] <- shared_id
  }
}

# Convert lists to matrices
if (length(lfc_list) == 0) stop("No comparisons extracted.")
lfc_matrix <- do.call(cbind, lfc_list)
se_matrix <- do.call(cbind, se_list)
pval_matrix <- do.call(cbind, pval_list)
rownames(lfc_matrix) <- rownames(cts)
rownames(se_matrix) <- rownames(cts)
rownames(pval_matrix) <- rownames(cts)
# Setup effective sample sizes (N_eff) vector needed for RE Stouffer weights
neff_vec <- unlist(neff_list)[colnames(pval_matrix)]

# --- 5. Metafor Heterogeneity Analysis ---
cat("\n=== RUNNING METAFOR HETEROGENEITY ANALYSIS ===\n")

cat("Computing tau^2, I^2, and meta-LFC (rma.mv) for", n_genes, "genes. This may take a moment...\n")

# Function to calculate Pseudo-I2 given sampling variances and tau2 component(s)
calculate_pseudo_i2 <- function(vi, tau2_vec) {
  W <- 1 / vi
  k <- length(vi)
  if (k > 1 && sum(W) > 0) {
    vt <- (k - 1) / (sum(W) - sum(W^2) / sum(W))
    total_var <- sum(tau2_vec, na.rm = TRUE) + vt
    return(100 * tau2_vec / total_var)
  } else {
    return(rep(NA, length(tau2_vec)))
  }
}

# Function to perform Cluster Wild Bootstrapping (manual equivalent of wildmeta)
run_manual_cwb <- function(orig_yi, orig_F, cluster_var, refit_fn, R = 399, null_preds = NULL, null_resids = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(null_preds) || is.null(null_resids)) {
    null_preds <- rep(0, length(orig_yi))
    null_resids <- orig_yi
  }

  unique_clusters <- unique(cluster_var)
  n_clusters <- length(unique_clusters)
  webb_dist <- c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5))

  boots <- numeric(R)
  for (b in seq_len(R)) {
    w <- sample(webb_dist, n_clusters, replace = TRUE)
    w_map <- setNames(w, unique_clusters)
    y_boot <- null_preds + null_resids * w_map[cluster_var]
    boots[b] <- refit_fn(y_boot)
  }

  valid_boots <- boots[!is.na(boots)]
  B_valid <- length(valid_boots)
  if (B_valid > 0) {
    p_val <- (sum(valid_boots > orig_F) + 1) / (B_valid + 1)
  } else {
    p_val <- NA
  }
  return(list(p_val = p_val, boots = valid_boots))
}


pb <- txtProgressBar(min = 0, max = n_genes, style = 3)

meta_results <- lapply(seq_len(n_genes), function(i) {
  setTxtProgressBar(pb, i)
  yi <- lfc_matrix[i, ]
  sei <- se_matrix[i, ]

  valid <- !is.na(yi) & !is.na(sei) & sei > 0 & !is.na(pval_matrix[i, ])
  valid_gse_all <- comparison_gse_map[valid]
  n_valid_datasets <- length(unique(valid_gse_all))

  if (n_valid_datasets >= max(2, min_ds)) {
    has_valid_results <- FALSE
    vi <- sei[valid]^2
    valid_gse <- valid_gse_all
    valid_comparison <- names(yi)[valid]

    # Fit Random-Effects Model (Multilevel)
    res_rma_re <- tryCatch(
      {
        metafor::rma.mv(yi = yi[valid], V = vi, random = ~ 1 | valid_gse, method = "REML", test = "t")
      },
      error = function(e) {
        tryCatch(
          {
            metafor::rma.mv(yi = yi[valid], V = vi, random = ~ 1 | valid_gse, method = "REML", test = "t", control = list(optimizer = "optim", optmethod = "Nelder-Mead"))
          },
          error = function(e2) NULL
        )
      }
    )

    # Exact variance-covariance matrix of dependent estimates
    V_mat <- tryCatch(
      {
        n_val <- length(valid_comparison)
        V_exact <- diag(vi, nrow = n_val, ncol = n_val)
        rownames(V_exact) <- valid_comparison
        colnames(V_exact) <- valid_comparison

        for (r in 1:n_val) {
          comp1 <- valid_comparison[r]
          gse1 <- valid_gse[r]
          for (c in 1:n_val) {
            comp2 <- valid_comparison[c]
            gse2 <- valid_gse[c]
            if (r != c && gse1 == gse2) {
              if (is.null(all_cov_mats[[gse1]])) {
                stop(sprintf("Missing covariance matrix for dataset '%s' during comparison '%s' vs '%s'", gse1, comp1, comp2))
              }

              true_cov <- NA
              if ("dataset_level" %in% names(all_cov_mats[[gse1]])) {
                true_cov <- all_cov_mats[[gse1]][["dataset_level"]][comp1, comp2]
              } else if (!is.null(all_cov_mats[[gse1]][[rownames(lfc_matrix)[i]]])) {
                true_cov <- all_cov_mats[[gse1]][[rownames(lfc_matrix)[i]]][comp1, comp2]
              }

              if (is.na(true_cov)) {
                stop(sprintf("Missing or NA covariance between '%s' and '%s' in dataset '%s' for gene '%s'", comp1, comp2, gse1, rownames(lfc_matrix)[i]))
              }

              sign1 <- sign_map[[comp1]]
              sign2 <- sign_map[[comp2]]
              if (is.null(sign1) || is.na(sign1) || is.null(sign2) || is.na(sign2)) {
                stop(sprintf("Missing contrast sign mapping for '%s' (sign=%s) or '%s' (sign=%s)", comp1, as.character(sign1), comp2, as.character(sign2)))
              }

              V_exact[r, c] <- true_cov * sign1 * sign2
            }
          }
        }
        V_exact <- (V_exact + t(V_exact)) / 2
        is_pd <- tryCatch(
          {
            chol(V_exact)
            TRUE
          },
          error = function(e) FALSE
        )
        if (!is_pd) {
          diag(V_exact) <- diag(V_exact) + 1e-6
        }
        V_exact
      },
      error = function(e) stop(paste("CRITICAL ERROR in V_mat for gene", rownames(lfc_matrix)[i], paste0("(index ", i, "):"), e$message))
    )

    # Fit 3-Level Random-Effects Model
    res_rma_3LRE <- tryCatch(
      {
        metafor::rma.mv(yi = yi[valid], V = V_mat, random = ~ 1 | valid_gse / valid_comparison, method = "REML")
      },
      error = function(e) {
        tryCatch(
          {
            metafor::rma.mv(yi = yi[valid], V = V_mat, random = ~ 1 | valid_gse / valid_comparison, method = "REML", control = list(optimizer = "optim", optmethod = "Nelder-Mead"))
          },
          error = function(e2) {
            if (i <= 10) cat("\n[Gene", i, "] 3LRE ERROR (all optimizers failed):", e2$message, "\n")
            NULL
          }
        )
      }
    )

    # Fit Fixed-Effects Model
    res_rma_fe <- tryCatch(
      {
        metafor::rma.mv(yi = yi[valid], V = V_mat, test = "t")
      },
      error = function(e) NULL
    )

    if (!is.null(res_rma_re)) {
      tau2 <- res_rma_re$sigma2[1]

      # Pseudo-I2 calculation
      I2 <- calculate_pseudo_i2(vi, tau2)

      re_lfc <- res_rma_re$beta[1]
      re_se <- res_rma_re$se[1]
      re_pval <- res_rma_re$pval[1]
      het_pval <- res_rma_re$QEp
    } else {
      re_lfc <- NA
      re_se <- NA
      re_pval <- NA
      tau2 <- NA
      I2 <- NA
      het_pval <- NA
    }

    if (!is.null(res_rma_3LRE)) {
      tau2_Level3 <- res_rma_3LRE$sigma2[1]
      tau2_Level2 <- res_rma_3LRE$sigma2[2]

      # Pseudo-I2 calculation
      i2_res <- calculate_pseudo_i2(vi, c(tau2_Level3, tau2_Level2))
      I2_Level3 <- i2_res[1]
      I2_Level2 <- i2_res[2]

      # Robust Variance Estimation (metafor + clubSandwich)
      robust_res <- tryCatch(
        {
          metafor::robust(res_rma_3LRE, cluster = valid_gse, clubSandwich = TRUE)
        },
        error = function(e) NULL
      )

      if (!is.null(robust_res)) {
        re3L_lfc <- as.numeric(robust_res$beta)[1]
        re3L_se <- as.numeric(robust_res$se)[1]
        re3L_pval <- as.numeric(robust_res$pval)[1]
        het3L_pval <- res_rma_3LRE$QEp

        # --- NEW: Cluster Wild Bootstrapping (Manual Implementation) ---
        boot_pval <- re3L_pval
        if (opt$use_cwb && !is.na(re3L_pval)) {
          # Fallback to the original p-value for genes that skip the bootstrap
          # so that Benjamini-Hochberg FDR correction doesn't artificially shrink
          boot_pval <- re3L_pval

          if (re3L_pval < opt$cwb_pval_thresh) {
            boot_pval <- tryCatch(
              {
                # Original robust Wald statistic (F-statistic for 1 parameter)
                orig_F <- as.numeric(robust_res$beta)[1]^2 / as.numeric(robust_res$se)[1]^2

                refit_fn <- function(y_boot) {
                  res_boot <- tryCatch(
                    {
                      metafor::rma.mv(
                        yi = y_boot, V = V_mat, random = ~ 1 | valid_gse / valid_comparison,
                        sigma2 = res_rma_3LRE$sigma2, method = "REML"
                      )
                    },
                    error = function(e) NULL
                  )

                  if (!is.null(res_boot)) {
                    rob_boot <- tryCatch(
                      {
                        metafor::robust(res_boot, cluster = valid_gse, clubSandwich = TRUE)
                      },
                      error = function(e) NULL
                    )
                    if (!is.null(rob_boot)) {
                      return(as.numeric(rob_boot$beta)[1]^2 / as.numeric(rob_boot$se)[1]^2)
                    }
                  }
                  return(NA)
                }

                cwb_res <- run_manual_cwb(
                  orig_yi = yi[valid],
                  orig_F = orig_F,
                  cluster_var = valid_gse,
                  refit_fn = refit_fn,
                  R = opt$cwb_r
                )

                cwb_res$p_val
              },
              error = function(e) NA
            )
          } # close if (re3L_pval < opt$cwb_pval_thresh)
        } # close if (opt$use_cwb && !is.na(re3L_pval))
      } else {
        # Fallback requested by user: NA if robust fails
        re3L_lfc <- NA
        re3L_se <- NA
        re3L_pval <- NA
        boot_pval <- NA
        tau2_Level3 <- NA
        tau2_Level2 <- NA
        I2_Level3 <- NA
        I2_Level2 <- NA
        het3L_pval <- NA
      }
    } else {
      re3L_lfc <- NA
      re3L_se <- NA
      re3L_pval <- NA
      boot_pval <- NA
      tau2_Level3 <- NA
      tau2_Level2 <- NA
      I2_Level3 <- NA
      I2_Level2 <- NA
      het3L_pval <- NA
    }

    if (!is.null(res_rma_fe)) {
      fe_lfc <- as.numeric(res_rma_fe$beta)[1]
      fe_se <- as.numeric(res_rma_fe$se)[1]
      fe_pval <- res_rma_fe$pval[1]
    } else {
      fe_lfc <- NA
      fe_se <- NA
      fe_pval <- NA
    }


    # --- NEW: Stouffer's Meta-Analysis (Two-Sided & Direction-Aware) ---
    stouffer_z <- NA
    stouffer_pval <- NA
    dir_stouffer_z <- NA
    dir_stouffer_pval <- NA
    re_stouffer_z <- NA
    re_stouffer_se <- NA
    re_stouffer_pval <- NA
    boot_re_stouffer_pval <- NA
    df_robust <- NA

    gene_pvals_stouffer <- pval_matrix[i, ]
    valid_stouffer <- !is.na(gene_pvals_stouffer)
    n_valid_datasets_stouffer <- length(unique(comparison_gse_map[valid_stouffer]))

    if (n_valid_datasets_stouffer >= max(2, min_ds)) {
      w_i <- sqrt(neff_vec[valid_stouffer])
      p_i <- gene_pvals_stouffer[valid_stouffer]

      n_stouffer <- sum(valid_stouffer)

      # Exact within-study correlation using precomputed Z-score covariance
      valid_cols <- colnames(pval_matrix)[valid_stouffer]
      valid_gse_stouffer <- comparison_gse_map[valid_cols]

      # Build working correlation matrix using exact correlations
      stouffer_mats <- tryCatch(
        {
          n_val <- length(valid_cols)
          V_exact_cor_dir <- diag(1, nrow = n_val, ncol = n_val)
          rownames(V_exact_cor_dir) <- valid_cols
          colnames(V_exact_cor_dir) <- valid_cols

          V_exact_cor_two <- diag(1, nrow = n_val, ncol = n_val)
          rownames(V_exact_cor_two) <- valid_cols
          colnames(V_exact_cor_two) <- valid_cols

          for (r in 1:n_val) {
            comp1 <- valid_cols[r]
            gse1 <- valid_gse_stouffer[r]
            for (c in 1:n_val) {
              comp2 <- valid_cols[c]
              gse2 <- valid_gse_stouffer[c]
              if (r != c && gse1 == gse2) {
                if (is.null(all_zcov_mats[[gse1]])) {
                  stop(sprintf("Missing Z-score covariance matrix for dataset '%s' during comparison '%s' vs '%s'", gse1, comp1, comp2))
                }

                cor_val <- NA
                if ("dataset_level" %in% names(all_zcov_mats[[gse1]])) {
                  cor_val <- all_zcov_mats[[gse1]][["dataset_level"]][comp1, comp2]
                } else if (!is.null(all_zcov_mats[[gse1]][[rownames(lfc_matrix)[i]]])) {
                  cor_val <- all_zcov_mats[[gse1]][[rownames(lfc_matrix)[i]]][comp1, comp2]
                }

                if (is.na(cor_val)) {
                  stop(sprintf("Missing or NA Z-score correlation between '%s' and '%s' in dataset '%s' for gene '%s'", comp1, comp2, gse1, rownames(lfc_matrix)[i]))
                }

                sign1 <- sign_map[[comp1]]
                sign2 <- sign_map[[comp2]]
                if (is.null(sign1) || is.na(sign1) || is.null(sign2) || is.na(sign2)) {
                  stop(sprintf("Missing contrast sign mapping for '%s' (sign=%s) or '%s' (sign=%s)", comp1, as.character(sign1), comp2, as.character(sign2)))
                }

                # Bound correlation to strictly within [-1 + 1e-6, 1 - 1e-6] to avoid floating-point issues
                cor_val <- max(-1 + 1e-6, min(1 - 1e-6, cor_val))

                V_exact_cor_two[r, c] <- cor_val
                V_exact_cor_dir[r, c] <- cor_val * sign1 * sign2
              }
            }
          }
          V_exact_cor_two <- (V_exact_cor_two + t(V_exact_cor_two)) / 2
          V_exact_cor_dir <- (V_exact_cor_dir + t(V_exact_cor_dir)) / 2

          is_pd_two <- tryCatch(
            {
              chol(V_exact_cor_two)
              TRUE
            },
            error = function(e) FALSE
          )
          if (!is_pd_two) diag(V_exact_cor_two) <- diag(V_exact_cor_two) + 1e-6

          is_pd_dir <- tryCatch(
            {
              chol(V_exact_cor_dir)
              TRUE
            },
            error = function(e) FALSE
          )
          if (!is_pd_dir) diag(V_exact_cor_dir) <- diag(V_exact_cor_dir) + 1e-6

          list(dir = V_exact_cor_dir, two = V_exact_cor_two)
        },
        error = function(e) stop(paste("CRITICAL ERROR in stouffer_mats for gene", rownames(lfc_matrix)[i], paste0("(index ", i, "):"), e$message))
      )

      V_stouffer_dir <- stouffer_mats$dir
      V_stouffer_twosided <- stouffer_mats$two

      # Variance of the weighted sum of correlated Z-scores: W^T * Cor * W
      var_sum_z_twosided <- as.numeric(t(w_i) %*% V_stouffer_twosided %*% w_i)
      var_sum_z_dir <- as.numeric(t(w_i) %*% V_stouffer_dir %*% w_i)

      # Avoid extreme p-values
      p_i[p_i < .Machine$double.eps] <- .Machine$double.eps
      p_i[p_i > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps

      # 1. Two-Sided Stouffer's (ignoring LFC direction)
      z_i_twosided <- qnorm(1 - p_i)
      if (n_stouffer > 1 && var_sum_z_twosided > 0) {
        stouffer_z <- sum(w_i * z_i_twosided) / sqrt(var_sum_z_twosided)
      } else {
        stouffer_z <- z_i_twosided[1]
      }
      stouffer_pval <- pnorm(stouffer_z, lower.tail = FALSE)

      # 2. Direction-Aware Stouffer's (incorporating LFC direction)
      lfc_i <- lfc_matrix[i, valid_stouffer]
      sign_i <- sign(lfc_i)
      sign_i[sign_i == 0] <- 1 # default to positive if LFC is exactly 0

      z_i_dir <- sign_i * qnorm(1 - p_i / 2)
      if (n_stouffer > 1 && var_sum_z_dir > 0) {
        dir_stouffer_z <- sum(w_i * z_i_dir) / sqrt(var_sum_z_dir)
      } else {
        dir_stouffer_z <- z_i_dir[1]
      }
      dir_stouffer_pval <- 2 * pnorm(-abs(dir_stouffer_z))

      # --- NEW: RE Stouffer's Meta-Analysis (rma.mv + robust) ---
      # We use Two-Sided Z-scores (z_i_twosided) as the dependent variable
      # Use assumed rho=0.5 working model for within-study correlation;
      # RVE (clubSandwich) corrects inference regardless of misspecification

      res_rma_re_stouffer <- tryCatch(
        {
          metafor::rma.mv(
            yi = z_i_twosided,
            V = V_stouffer_twosided,
            W = neff_vec[valid_stouffer],
            random = ~ 1 | valid_gse_stouffer / valid_cols,
            method = "REML"
          )
        },
        error = function(e) {
          # Retry with Nelder-Mead optimizer if the default nlminb fails to converge
          tryCatch(
            {
              metafor::rma.mv(
                yi = z_i_twosided,
                V = V_stouffer_twosided,
                W = neff_vec[valid_stouffer],
                random = ~ 1 | valid_gse_stouffer / valid_cols,
                method = "REML",
                control = list(optimizer = "optim", optmethod = "Nelder-Mead")
              )
            },
            error = function(e2) {
              if (i <= 10) cat("\n[Gene", i, "] REStouffer ERROR (all optimizers failed):", e2$message, "\n")
              NULL
            }
          )
        }
      )

      if (!is.null(res_rma_re_stouffer)) {
        robust_re_stouffer <- tryCatch(
          {
            metafor::robust(res_rma_re_stouffer, cluster = valid_gse_stouffer, clubSandwich = TRUE)
          },
          error = function(e) NULL
        )

        if (!is.null(robust_re_stouffer)) {
          re_stouffer_z <- as.numeric(robust_re_stouffer$beta)[1]
          re_stouffer_se <- as.numeric(robust_re_stouffer$se)[1]

          # One-tailed test: testing if the pooled Z-score > 0
          t_stat <- re_stouffer_z / re_stouffer_se
          df_robust <- robust_re_stouffer$ddf[1]
          if (i <= 5) {
            cat(
              "\n[Gene", i, "] Full ddf array:", paste(robust_re_stouffer$ddf, collapse = ", "),
              " | Extracted df_robust:", df_robust, "\n"
            )
          }
          if (!is.na(df_robust) && df_robust > 0) {
            re_stouffer_pval <- pt(t_stat, df = df_robust, lower.tail = FALSE)
          } else {
            re_stouffer_pval <- NA
          }

          # --- NEW: Cluster Wild Bootstrapping for RE Stouffer ---
          boot_re_stouffer_pval <- re_stouffer_pval
          if (opt$use_cwb_stouffer && !is.na(re_stouffer_pval)) {
            if (re_stouffer_pval < opt$cwb_stouffer_pval_thresh) {
              boot_re_stouffer_pval <- tryCatch(
                {
                  orig_F_stouffer <- t_stat^2
                  refit_fn_stouffer <- function(y_boot) {
                    res_b <- tryCatch(
                      {
                        metafor::rma.mv(
                          yi = y_boot, V = V_stouffer, W = neff_vec[valid_stouffer],
                          random = ~ 1 | valid_gse_stouffer / valid_cols, method = "REML"
                        )
                      },
                      error = function(e) NULL
                    )
                    if (!is.null(res_b)) {
                      rob_b <- tryCatch(metafor::robust(res_b, cluster = valid_gse_stouffer, clubSandwich = TRUE), error = function(e) NULL)
                      if (!is.null(rob_b)) {
                        return(as.numeric(rob_b$beta)[1]^2 / as.numeric(rob_b$se)[1]^2)
                      }
                    }
                    return(NA)
                  }
                  cwb_res_stouffer <- run_manual_cwb(
                    orig_yi = z_i_twosided, orig_F = orig_F_stouffer,
                    cluster_var = valid_gse_stouffer, refit_fn = refit_fn_stouffer,
                    R = opt$cwb_stouffer_r
                  )
                  if (!is.na(cwb_res_stouffer$p_val)) {
                    if (re_stouffer_z > 0) {
                      cwb_res_stouffer$p_val / 2
                    } else {
                      1 - (cwb_res_stouffer$p_val / 2)
                    }
                  } else {
                    NA
                  }
                },
                error = function(e) NA
              )
            }
          }
        }
      }
    }



    if (has_valid_results || !is.null(res_rma_re) || !is.null(res_rma_fe) || !is.null(res_rma_3LRE)) {
      return(c(
        Meta_3LRE_LFC = re3L_lfc, Meta_3LRE_SE = re3L_se, Meta_3LRE_pval = re3L_pval,
        Meta_3LRE_boot_pval = boot_pval,
        Meta_FE_LFC = fe_lfc, Meta_FE_SE = fe_se, Meta_FE_pval = fe_pval,
        Meta_Stouffer_Z = stouffer_z,
        Meta_Stouffer_pval = stouffer_pval,
        Meta_DirectionalStouffer_Z = dir_stouffer_z,
        Meta_DirectionalStouffer_pval = dir_stouffer_pval,
        Meta_RE_Stouffer_Z = re_stouffer_z, Meta_RE_Stouffer_SE = re_stouffer_se, Meta_RE_Stouffer_pval = re_stouffer_pval, Meta_RE_Stouffer_boot_pval = boot_re_stouffer_pval, Meta_RE_Stouffer_ddf = df_robust,
        Heterogeneity_tau2 = tau2, Heterogeneity_I2 = I2, Heterogeneity_pval = het_pval,
        Heterogeneity_3L_tau2_Level3 = tau2_Level3, Heterogeneity_3L_tau2_Level2 = tau2_Level2,
        Heterogeneity_3L_I2_Level3 = I2_Level3, Heterogeneity_3L_I2_Level2 = I2_Level2,
        Heterogeneity_3L_pval = het3L_pval
      ))
    }
  }
  return(c(
    Meta_3LRE_LFC = NA, Meta_3LRE_SE = NA, Meta_3LRE_pval = NA,
    Meta_3LRE_boot_pval = NA,
    Meta_FE_LFC = NA, Meta_FE_SE = NA, Meta_FE_pval = NA,
    Meta_Stouffer_Z = NA,
    Meta_Stouffer_pval = NA,
    Meta_DirectionalStouffer_Z = NA,
    Meta_DirectionalStouffer_pval = NA,
    Meta_RE_Stouffer_Z = NA, Meta_RE_Stouffer_SE = NA, Meta_RE_Stouffer_pval = NA, Meta_RE_Stouffer_boot_pval = NA, Meta_RE_Stouffer_ddf = NA,
    Heterogeneity_tau2 = NA, Heterogeneity_I2 = NA, Heterogeneity_pval = NA,
    Heterogeneity_3L_tau2_Level3 = NA, Heterogeneity_3L_tau2_Level2 = NA,
    Heterogeneity_3L_I2_Level3 = NA, Heterogeneity_3L_I2_Level2 = NA,
    Heterogeneity_3L_pval = NA
  ))
})
close(pb)

meta_results_df <- as.data.frame(do.call(rbind, meta_results))
meta_results_df$cluster_id <- rownames(lfc_matrix)
meta_results_df$Meta_3LRE_padj <- p.adjust(meta_results_df$Meta_3LRE_pval, method = "BH")
meta_results_df$Meta_3LRE_boot_padj <- p.adjust(meta_results_df$Meta_3LRE_boot_pval, method = "BH")
meta_results_df$Meta_FE_padj <- p.adjust(meta_results_df$Meta_FE_pval, method = "BH")
meta_results_df$Meta_Stouffer_padj <- p.adjust(meta_results_df$Meta_Stouffer_pval, method = "BH")
meta_results_df$Meta_DirectionalStouffer_padj <- p.adjust(meta_results_df$Meta_DirectionalStouffer_pval, method = "BH")
meta_results_df$Meta_RE_Stouffer_padj <- p.adjust(meta_results_df$Meta_RE_Stouffer_pval, method = "BH")
meta_results_df$Meta_RE_Stouffer_boot_padj <- p.adjust(meta_results_df$Meta_RE_Stouffer_boot_pval, method = "BH")
meta_results_df$Heterogeneity_padj <- p.adjust(meta_results_df$Heterogeneity_pval, method = "BH")
meta_results_df$Heterogeneity_3L_padj <- p.adjust(meta_results_df$Heterogeneity_3L_pval, method = "BH")

# Variance-weighted mean of absolute LFCs (Inverse Variance Weighting)
W_mat <- 1 / (se_matrix^2)
W_mat[is.na(lfc_matrix) | is.na(se_matrix) | se_matrix <= 0] <- NA
meta_results_df$Mean_Absolute_LFC <- rowSums(W_mat * abs(lfc_matrix), na.rm = TRUE) / rowSums(W_mat, na.rm = TRUE)

# Calculate Average BaseMean across datasets for each gene (for MA plot)
baseMean_matrix <- matrix(NA, nrow = n_genes, ncol = length(unique_gses))
rownames(baseMean_matrix) <- rownames(lfc_matrix)
colnames(baseMean_matrix) <- unique_gses

for (gse_id in unique_gses) {
  if (!is.null(dataset_dds_list[[gse_id]])) {
    res_tmp <- results(dataset_dds_list[[gse_id]])
    baseMean_matrix[rownames(res_tmp), gse_id] <- res_tmp$baseMean
  }
}
avg_baseMean <- rowMeans(baseMean_matrix, na.rm = TRUE)
meta_results_df$baseMean <- avg_baseMean


# Export Results
res_df_export <- data.frame(cluster_id = rownames(lfc_matrix), stringsAsFactors = FALSE)
res_df_export <- left_join(res_df_export, cluster_annotation, by = "cluster_id")
res_df_export <- left_join(res_df_export, meta_results_df, by = "cluster_id")

res_df_export <- res_df_export %>% select(
  cluster_id, Representative_Gene, All_Genes_Locus_Tags, baseMean, Mean_Absolute_LFC,
  Meta_Stouffer_Z, Meta_Stouffer_pval, Meta_Stouffer_padj,
  Meta_DirectionalStouffer_Z, Meta_DirectionalStouffer_pval, Meta_DirectionalStouffer_padj,
  Meta_RE_Stouffer_Z, Meta_RE_Stouffer_SE, Meta_RE_Stouffer_pval, Meta_RE_Stouffer_padj,
  Meta_3LRE_LFC, Meta_3LRE_SE, Meta_3LRE_pval, Meta_3LRE_padj,
  Meta_3LRE_boot_pval, Meta_3LRE_boot_padj,
  Meta_FE_LFC, Meta_FE_SE, Meta_FE_pval, Meta_FE_padj,
  Heterogeneity_3L_I2_Level3, Heterogeneity_3L_I2_Level2, Heterogeneity_3L_tau2_Level3, Heterogeneity_3L_tau2_Level2, Heterogeneity_3L_pval, Heterogeneity_3L_padj,
  Heterogeneity_I2, Heterogeneity_tau2, Heterogeneity_pval, Heterogeneity_padj,
  everything()
)

n_analyzed <- sum(!is.na(res_df_export$Meta_RE_Stouffer_pval))
n_total_genes <- nrow(res_df_export)
cat("\n=== META-ANALYSIS SUMMARY ===\n")
cat(sprintf(
  "  %d out of %d genes (%.1f%%) had enough datasets (>= %d) to be meta-analyzed.\n\n",
  n_analyzed, n_total_genes, 100 * n_analyzed / n_total_genes, max(2, min_ds)
))


out_file_all <- file.path(out_dir, "DESeq2_contrast_results.csv")
write.csv(res_df_export, file = out_file_all, row.names = FALSE)
cat("✓ Exported", out_file_all, "\n")

# Export Significant Results - 3-Level RE Model
res_df_sig_3lre <- res_df_export %>% filter(!is.na(Meta_3LRE_padj) & Meta_3LRE_padj < 0.05 & abs(Meta_3LRE_LFC) > 0.5)
write.csv(res_df_sig_3lre, file = file.path(out_dir, "DESeq2_contrast_results_significant_Meta3LRE.csv"), row.names = FALSE)

# Export Significant Results - 3-Level RE Model (Bootstrapped)
if (opt$use_cwb) {
  res_df_sig_3lre_boot <- res_df_export %>% filter(!is.na(Meta_3LRE_boot_padj) & Meta_3LRE_boot_padj < 0.05 & abs(Meta_3LRE_LFC) > 0.5)
  write.csv(res_df_sig_3lre_boot, file = file.path(out_dir, "DESeq2_contrast_results_significant_Meta3LRE_Boot.csv"), row.names = FALSE)
}

# Export Significant Results - FE Model
res_df_sig_fe <- res_df_export %>% filter(!is.na(Meta_FE_padj) & Meta_FE_padj < 0.05 & abs(Meta_FE_LFC) > 0.5)
write.csv(res_df_sig_fe, file = file.path(out_dir, "DESeq2_contrast_results_significant_MetaFE.csv"), row.names = FALSE)

# Export Significant Results - Two-Sided Stouffer's
res_df_sig_stouffer <- res_df_export %>% filter(!is.na(Meta_Stouffer_padj) & Meta_Stouffer_padj < 0.05 & !is.na(Mean_Absolute_LFC) & Mean_Absolute_LFC > 0.5)
write.csv(res_df_sig_stouffer, file = file.path(out_dir, "DESeq2_contrast_results_significant_MetaStouffer.csv"), row.names = FALSE)

# Export Significant Results - Direction-Aware Stouffer's
res_df_sig_dir_stouffer <- res_df_export %>% filter(!is.na(Meta_DirectionalStouffer_padj) & Meta_DirectionalStouffer_padj < 0.05 & !is.na(Meta_FE_LFC) & Meta_FE_LFC > 0.5)
write.csv(res_df_sig_dir_stouffer, file = file.path(out_dir, "DESeq2_contrast_results_significant_MetaDirectionalStouffer.csv"), row.names = FALSE)

# Export Significant Results - RE Stouffer's
res_df_sig_re_stouffer <- res_df_export %>% filter(!is.na(Meta_RE_Stouffer_padj) & Meta_RE_Stouffer_padj < 0.05 & !is.na(Mean_Absolute_LFC) & Mean_Absolute_LFC > 0.5)
write.csv(res_df_sig_re_stouffer, file = file.path(out_dir, "DESeq2_contrast_results_significant_MetaREStouffer.csv"), row.names = FALSE)



# Export Individual Dataset Coefficients (Sign-Adjusted LFCs and SEs)
coef_df_export <- data.frame(cluster_id = rownames(lfc_matrix), stringsAsFactors = FALSE)
coef_df_export <- left_join(coef_df_export, cluster_annotation, by = "cluster_id")

for (comp_name in colnames(lfc_matrix)) {
  lfc_col <- paste0(comp_name, "_LFC_adj")
  se_col <- paste0(comp_name, "_SE")

  coef_df_export[[lfc_col]] <- lfc_matrix[, comp_name]
  coef_df_export[[se_col]] <- se_matrix[, comp_name]
}

out_file_coefs <- file.path(out_dir, "dataset_coefficients.csv")
write.csv(coef_df_export, file = out_file_coefs, row.names = FALSE)
cat("✓ Exported dataset coefficients to", out_file_coefs, "\n")

# --- TopGO / ORA Meta-Analysis ---
source("ora_meta_analysis.R")
run_ora_meta_analysis(
  gaf_file = opt$gaf_file,
  lfc_list = lfc_list,
  pval_list = pval_list,
  cts = cts,
  cluster_annotation = cluster_annotation,
  comparison_gse_map = comparison_gse_map,
  min_ds = min_ds,
  out_dir = out_dir,
  ora_padj = opt$ora_padj,
  ora_lfc = opt$ora_lfc
)

# --- 6. Visualizations ---
cat("\n=== GENERATING VISUALIZATIONS ===\n")

# Volcano Plot Generator Function
generate_volcano_plot <- function(res_df_export, padj_col, lfc_col, title, x_label, y_label, output_file, vis_dir) {
  cat(paste0("Generating Volcano Plot: ", title, "...\n"))
  volc_df <- res_df_export
  volc_df$Significant <- ifelse(!is.na(volc_df[[padj_col]]) & volc_df[[padj_col]] < 0.05 & abs(volc_df[[lfc_col]]) > 0.5, "Significant", "Not Significant")

  top_genes <- volc_df %>%
    filter(Significant == "Significant") %>%
    arrange(.data[[padj_col]]) %>%
    head(20)

  volc_df$Label <- ifelse(volc_df$cluster_id %in% top_genes$cluster_id, volc_df$Representative_Gene, "")

  p_volcano <- ggplot(volc_df, aes(x = .data[[lfc_col]], y = -log10(.data[[padj_col]]), color = Significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_text_repel(aes(label = Label), size = 3.5, max.overlaps = 50, show.legend = FALSE, color = "black") +
    theme_minimal() +
    labs(
      title = title,
      x = x_label, y = y_label
    )

  ggsave(file.path(vis_dir, output_file), p_volcano, width = 8, height = 6, dpi = 300)
  cat(paste0("✓ ", output_file, " saved\n"))
}

# Generate Volcano Plots for each model
generate_volcano_plot(res_df_export, "Meta_FE_padj", "Meta_FE_LFC", "Volcano Plot - Fixed Effects Meta-Analysis", "Fixed Effects Log2 Fold Change", "-Log10(FE p-adj)", "Volcano_Plot_FE.png", vis_dir)
generate_volcano_plot(res_df_export, "Meta_3LRE_padj", "Meta_3LRE_LFC", "Volcano Plot - 3-Level Random Effects Meta-Analysis", "3-Level Random Effects Log2 Fold Change", "-Log10(3L-RE p-adj)", "Volcano_Plot_3LRE.png", vis_dir)
if (opt$use_cwb) {
  generate_volcano_plot(res_df_export, "Meta_3LRE_boot_padj", "Meta_3LRE_LFC", "Volcano Plot - 3-Level RE with Cluster Wild Bootstrap", "3-Level Random Effects Log2 Fold Change", "-Log10(3L-RE Bootstrap p-adj)", "Volcano_Plot_3LRE_Boot.png", vis_dir)
}


# P-Value Histogram (Fixed Effects)
cat("Creating p-value histogram (Fixed Effects)...\n")
tryCatch(
  {
    png(file.path(vis_dir, "PValue_Histogram_FE.png"), width = 800, height = 600, res = 150)
    hist(res_df_export$Meta_FE_pval,
      breaks = 50, col = "grey",
      main = "Fixed Effects p-Value Distribution",
      xlab = "Unadjusted FE p-value", ylab = "Frequency"
    )
    dev.off()
    cat("✓ p-Value histogram saved\n")
  },
  error = function(e) {
    cat("⚠ Warning: Could not create p-value histogram:", e$message, "\n")
  }
)

# MA Plot (Fixed Effects)
cat("Creating MA plot (Fixed Effects)...\n")
tryCatch(
  {
    png(file.path(vis_dir, "MA_plot_FE.png"), width = 8, height = 6, units = "in", res = 300)
    plot_data <- res_df_export[!is.na(res_df_export$baseMean) & !is.na(res_df_export$Meta_FE_LFC), ]

    plot_data$col <- ifelse(!is.na(plot_data$Meta_FE_padj) & plot_data$Meta_FE_padj < 0.05, "blue", "grey")

    y_max <- min(3, max(abs(plot_data$Meta_FE_LFC), na.rm = TRUE))
    plot_data$Meta_FE_LFC_capped <- pmax(pmin(plot_data$Meta_FE_LFC, y_max), -y_max)

    plot(plot_data$baseMean, plot_data$Meta_FE_LFC_capped,
      log = "x", pch = 16, cex = 0.5, col = plot_data$col,
      xlab = "mean of normalized counts (average across datasets)",
      ylab = "Fixed Effects Log2 Fold Change",
      main = "MA Plot - Fixed Effects Meta-Analysis",
      ylim = c(-y_max, y_max)
    )
    abline(h = 0, col = "red", lwd = 2)
    dev.off()
    cat("✓ MA plot saved\n")
  },
  error = function(e) {
    cat("⚠ Warning: Could not create MA plot:", e$message, "\n")
  }
)


# Function to generate Forest Plots
generate_forest_plot <- function(res_df_export, lfc_matrix, se_matrix, model_type, vis_dir) {
  if (model_type == "3LRE") {
    padj_col <- "Meta_3LRE_padj"
    lfc_col <- "Meta_3LRE_LFC"
    se_col <- "Meta_3LRE_SE"
    overall_label <- "Overall (3-Level RE)"
    plot_title <- "Forest Plot of Log2 Fold Changes (Top 15 3-Level RE Significant Genes)"
    file_name <- "ForestPlot_Top15_3LRE_Genes.png"
  } else if (model_type == "3LRE_Boot") {
    padj_col <- "Meta_3LRE_boot_padj"
    lfc_col <- "Meta_3LRE_LFC"
    se_col <- "Meta_3LRE_SE"
    overall_label <- "Overall (3-Level RE Bootstrap)"
    plot_title <- "Forest Plot of Log2 Fold Changes (Top 15 3-Level RE Bootstrapped Genes)"
    file_name <- "ForestPlot_Top15_3LRE_Boot_Genes.png"
  } else if (model_type == "FE") {
    padj_col <- "Meta_FE_padj"
    lfc_col <- "Meta_FE_LFC"
    se_col <- "Meta_FE_SE"
    overall_label <- "Overall (Fixed Effects)"
    plot_title <- "Forest Plot of Log2 Fold Changes (Top 15 FE Significant Genes)"
    file_name <- "ForestPlot_Top15_FE_Genes.png"
  } else if (model_type == "Stouffer") {
    padj_col <- "Meta_Stouffer_padj"
    lfc_col <- "Mean_Absolute_LFC"
    se_col <- NA
    overall_label <- "Overall (Stouffer)"
    plot_title <- "Forest Plot of Log2 Fold Changes (Top 15 Stouffer Significant Genes)"
    file_name <- "ForestPlot_Top15_Stouffer_Genes.png"
  } else if (model_type == "DirectionalStouffer") {
    padj_col <- "Meta_DirectionalStouffer_padj"
    lfc_col <- "Meta_FE_LFC"
    se_col <- NA
    overall_label <- "Overall (Dir. Stouffer)"
    plot_title <- "Forest Plot of Log2 Fold Changes (Top 15 Dir. Stouffer Significant Genes)"
    file_name <- "ForestPlot_Top15_DirStouffer_Genes.png"
  } else if (model_type == "RE_Stouffer") {
    padj_col <- "Meta_RE_Stouffer_padj"
    lfc_col <- "Mean_Absolute_LFC"
    se_col <- NA
    overall_label <- "Overall (RE Stouffer)"
    plot_title <- "Forest Plot of Log2 Fold Changes (Top 15 RE Stouffer Significant Genes)"
    file_name <- "ForestPlot_Top15_REStouffer_Genes.png"
  } else {
    stop("Unknown model type")
  }

  cat(paste0("Generating Forest Plot of LFCs for top 15 ", model_type, " genes...\n"))
  tryCatch(
    {
      top15 <- res_df_export %>%
        filter(!is.na(.data[[padj_col]]) & .data[[padj_col]] < 0.05 & abs(.data[[lfc_col]]) > 0.5) %>%
        arrange(.data[[padj_col]]) %>%
        head(15)

      if (nrow(top15) > 0) {
        top15_ids <- top15$cluster_id

        plot_data_list <- list()
        for (gene_id in top15_ids) {
          gene_name <- top15$Representative_Gene[top15$cluster_id == gene_id]
          gene_label <- paste0(gene_name, " (", gene_id, ")")

          # Individual comparisons
          comps <- colnames(lfc_matrix)
          gene_lfcs <- lfc_matrix[gene_id, comps]
          gene_ses <- se_matrix[gene_id, comps]

          valid <- !is.na(gene_lfcs) & !is.na(gene_ses)
          comps_valid <- comps[valid]

          if (length(comps_valid) > 0) {
            df_ind <- data.frame(
              Gene = gene_label,
              Comparison = comps_valid,
              Type = "Individual",
              LFC = gene_lfcs[valid],
              SE = gene_ses[valid],
              stringsAsFactors = FALSE
            )
            plot_data_list[[paste0(gene_id, "_ind")]] <- df_ind
          }

          # Meta-analysis summary point
          meta_lfc <- if (!is.na(lfc_col) && lfc_col %in% colnames(top15)) top15[[lfc_col]][top15$cluster_id == gene_id] else NA
          meta_se <- if (!is.na(se_col) && se_col %in% colnames(top15)) top15[[se_col]][top15$cluster_id == gene_id] else NA

          if (!is.null(meta_lfc) && !is.na(meta_lfc) && !is.null(meta_se) && !is.na(meta_se)) {
            df_meta <- data.frame(
              Gene = gene_label,
              Comparison = overall_label,
              Type = "Meta-Analysis",
              LFC = meta_lfc,
              SE = meta_se,
              stringsAsFactors = FALSE
            )
            plot_data_list[[paste0(gene_id, "_meta")]] <- df_meta
          }
        }

        plot_data <- do.call(rbind, plot_data_list)
        if (!is.null(plot_data) && nrow(plot_data) > 0) {
          comps_unique <- setdiff(unique(plot_data$Comparison), overall_label)
          plot_data$Comparison <- factor(plot_data$Comparison, levels = c(rev(comps_unique), overall_label))

          ordered_genes <- unique(paste0(top15$Representative_Gene, " (", top15$cluster_id, ")"))
          plot_data$Gene <- factor(plot_data$Gene, levels = ordered_genes)

          p_forest <- ggplot(plot_data, aes(x = LFC, y = Comparison, color = Type)) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
            geom_point(aes(shape = Type, size = Type)) +
            geom_errorbarh(aes(xmin = LFC - 1.96 * SE, xmax = LFC + 1.96 * SE), height = 0.2) +
            scale_color_manual(values = c("Individual" = "black", "Meta-Analysis" = "red")) +
            scale_shape_manual(values = c("Individual" = 16, "Meta-Analysis" = 18)) +
            scale_size_manual(values = c("Individual" = 2, "Meta-Analysis" = 3)) +
            facet_wrap(~Gene, scales = "free_y", ncol = 3) +
            theme_bw() +
            theme(
              legend.position = "bottom",
              strip.text = element_text(face = "bold")
            ) +
            labs(
              title = plot_title,
              x = "Log2 Fold Change (± 95% CI)", y = ""
            )

          ggsave(file.path(vis_dir, file_name), p_forest, width = 14, height = 12, dpi = 300)
          cat(paste0("✓ Forest Plot (", model_type, ") saved\n"))
        } else {
          cat(paste0("⚠ No data to plot for Forest Plot (", model_type, ").\n"))
        }
      } else {
        cat(paste0("⚠ No significant ", model_type, " genes to plot.\n"))
      }
    },
    error = function(e) {
      cat(paste0("⚠ Warning: Could not create Forest Plot (", model_type, "): "), e$message, "\n")
    }
  )
}

# Generate plots
generate_forest_plot(res_df_export, lfc_matrix, se_matrix, "3LRE", vis_dir)
if (opt$use_cwb) {
  generate_forest_plot(res_df_export, lfc_matrix, se_matrix, "3LRE_Boot", vis_dir)
}
generate_forest_plot(res_df_export, lfc_matrix, se_matrix, "FE", vis_dir)
generate_forest_plot(res_df_export, lfc_matrix, se_matrix, "Stouffer", vis_dir)
generate_forest_plot(res_df_export, lfc_matrix, se_matrix, "DirectionalStouffer", vis_dir)
generate_forest_plot(res_df_export, lfc_matrix, se_matrix, "RE_Stouffer", vis_dir)

# Profile Likelihood Plot for Top 3LRE Gene
cat("\nGenerating Profile Likelihood Plot for the top 3-Level RE gene...\n")
tryCatch(
  {
    top_3lre_gene <- res_df_export %>%
      filter(!is.na(Meta_3LRE_padj)) %>%
      arrange(Meta_3LRE_padj) %>%
      head(1)

    if (nrow(top_3lre_gene) > 0) {
      gene_id <- top_3lre_gene$cluster_id
      gene_name <- top_3lre_gene$Representative_Gene

      yi <- lfc_matrix[gene_id, ]
      sei <- se_matrix[gene_id, ]
      valid <- !is.na(yi) & !is.na(sei) & sei > 0
      valid_gse_all <- comparison_gse_map[valid]
      n_valid_datasets <- length(unique(valid_gse_all))

      if (n_valid_datasets >= max(2, min_ds)) {
        vi <- sei[valid]^2
        valid_gse <- comparison_gse_map[valid]
        valid_comparison <- names(yi)[valid]

        V_mat_top <- tryCatch(
          {
            metafor::vcalc(vi = vi, cluster = valid_gse, obs = valid_comparison, rho = 0.5)
          },
          error = function(e) vi
        )

        res_rma_3LRE_top <- tryCatch(
          {
            metafor::rma.mv(yi = yi[valid], V = V_mat_top, random = ~ 1 | valid_gse / valid_comparison, method = "REML")
          },
          error = function(e) NULL
        )

        if (!is.null(res_rma_3LRE_top)) {
          png(file.path(vis_dir, paste0("Profile_Likelihood_Top_3LRE_Gene.png")), width = 800, height = 400, res = 100)
          par(mfrow = c(1, 2))
          profile(res_rma_3LRE_top, main = paste("Profile Likelihood:", gene_name))
          dev.off()
          cat("✓ Profile Likelihood plot saved for top gene:", gene_name, "(", gene_id, ")\n")
        }
      }
    }
  },
  error = function(e) {
    cat("⚠ Warning: Could not create Profile Likelihood plot:", e$message, "\n")
  }
)

# --- NEW: Plot DDF Distribution ---
ddf_values <- meta_results_df$Meta_RE_Stouffer_ddf[!is.na(meta_results_df$Meta_RE_Stouffer_ddf)]
if (length(ddf_values) > 0) {
  tryCatch(
    {
      png(file.path(vis_dir, "RE_Stouffer_DDF_Distribution.png"), width = 800, height = 600, res = 150)
      boxplot(ddf_values,
        col = "skyblue", border = "darkblue",
        horizontal = TRUE,
        main = "Distribution of Satterthwaite ddf\n(RE Stouffer Model)",
        xlab = "Degrees of Freedom (ddf)"
      )
      dev.off()
      cat("\n✓ DDF Distribution plot saved to:", file.path(vis_dir, "RE_Stouffer_DDF_Distribution.png"), "\n")
    },
    error = function(e) {
      cat("⚠ Warning: Could not create DDF distribution plot:", e$message, "\n")
    }
  )
}

# --- NEW: Plot Z-score Input Distributions ---
tryCatch(
  {
    cat("\n=== PLOTTING Z-SCORE INPUT DISTRIBUTIONS ===\n")

    # Flatten matrices to get all individual inputs
    all_p <- as.numeric(pval_matrix)
    all_lfc <- as.numeric(lfc_matrix)

    valid_idx <- !is.na(all_p) & !is.na(all_lfc)
    p_valid <- all_p[valid_idx]
    lfc_valid <- all_lfc[valid_idx]

    # Avoid extreme p-values matching the loop logic
    p_valid[p_valid < .Machine$double.eps] <- .Machine$double.eps
    p_valid[p_valid > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps

    # Re-calculate individual z_i components
    z_twosided <- qnorm(1 - p_valid)

    sign_valid <- sign(lfc_valid)
    sign_valid[sign_valid == 0] <- 1
    z_dir <- sign_valid * qnorm(1 - p_valid / 2)

    # Function to plot histogram with overlay
    plot_z_dist <- function(z_vals, title, filename) {
      png(file.path(vis_dir, filename), width = 800, height = 600, res = 150)

      # Remove infs if any
      z_vals <- z_vals[is.finite(z_vals)]

      hist(z_vals,
        breaks = 100, freq = FALSE,
        col = "lightblue", border = "white",
        main = title, xlab = "Z-score"
      )

      # Overlay theoretical N(0,1)
      x_range <- seq(min(z_vals), max(z_vals), length.out = 500)
      lines(x_range, dnorm(x_range, mean = 0, sd = 1), col = "red", lwd = 2)

      # Overlay empirical density
      lines(density(z_vals), col = "blue", lwd = 2, lty = 2)

      legend("topright",
        legend = c("Theoretical N(0,1)", "Empirical Density"),
        col = c("red", "blue"), lwd = 2, lty = c(1, 2)
      )

      dev.off()
    }

    plot_z_dist(z_twosided, "Distribution of z_i_twosided inputs", "Z_Score_Dist_TwoSided.png")
    plot_z_dist(z_dir, "Distribution of z_i_dir inputs", "Z_Score_Dist_Directional.png")

    cat("✓ Z-score distribution plots saved to:", vis_dir, "\n")
  },
  error = function(e) {
    cat("⚠ Warning: Could not create Z-score distribution plots:", e$message, "\n")
  }
)


cat("\n=== ANALYSIS COMPLETE ===\n")
sink()
