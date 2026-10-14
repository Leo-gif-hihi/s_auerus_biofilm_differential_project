# dropout_test.R - Drop-Out / Sensitivity Permutation Analysis for Meta-Analysis
#
# This module performs leave-one-out and combinatorial dataset drop-out tests
# across all loaded GSE datasets (RNA-seq and microarray) down to a specified minimum.
# For each dataset combination, it runs 3-Level Random Effects (3LRE with clubSandwich RVE)
# and Fixed Effects (FE) meta-analysis models, outputs detailed CSVs per combination,
# and compiles a master dropout_summary.csv.

suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
  library(clubSandwich)
})

# Helper function to calculate Pseudo-I2 given sampling variances and tau2 component(s)
calculate_pseudo_i2_dropout <- function(vi, tau2_vec) {
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

run_dropout_test <- function(lfc_matrix,
                              se_matrix,
                              pval_matrix,
                              comparison_gse_map,
                              sign_map,
                              all_cov_mats,
                              cluster_annotation,
                              min_datasets = 0,
                              min_dropout_datasets = 3,
                              out_dir = ".",
                              force_rerun = FALSE) {

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("=== RUNNING DATASET DROP-OUT / SENSITIVITY TEST ===\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")

  dropout_dir <- file.path(out_dir, "drop_out_test")
  if (!dir.exists(dropout_dir)) {
    dir.create(dropout_dir, recursive = TRUE)
  }

  all_gses <- sort(unique(as.character(comparison_gse_map)))
  all_gses <- all_gses[!is.na(all_gses) & all_gses != ""]
  N <- length(all_gses)

  min_k <- max(2, as.integer(min_dropout_datasets))

  cat(sprintf("Total unique GSE datasets available: %d (%s)\n", N, paste(all_gses, collapse = ", ")))
  cat(sprintf("Configured minimum GSE datasets per combination (--min_dropout_datasets): %d\n", min_k))
  if (force_rerun) cat("  ⚠ --force_rerun_dropout enabled: Existing combination results will be recalculated.\n")

  if (N < min_k) {
    cat(sprintf("  ⚠ Warning: Total GSE datasets (%d) is less than min_dropout_datasets (%d). Cannot perform drop-out test.\n", N, min_k))
    return(invisible(NULL))
  }

  if (N <= 2) {
    cat("  ⚠ Notice: With 2 or fewer total datasets, combinatorial drop-out is not applicable.\n")
    return(invisible(NULL))
  }

  # Generate all combinations of size k from N-1 down to min_k
  comb_list <- list()
  for (k in (N - 1):min_k) {
    k_combs <- utils::combn(all_gses, k, simplify = FALSE)
    comb_list <- c(comb_list, k_combs)
  }

  n_combs <- length(comb_list)
  cat(sprintf("Total combinations to evaluate: %d\n\n", n_combs))

  all_genes <- rownames(lfc_matrix)
  n_genes <- length(all_genes)

  # Internal function to fit 3LRE & FE for any subset of GSE datasets
  fit_one_combination <- function(comb_gses, comb_name, comb_dir = NULL) {
    k_cur <- length(comb_gses)
    dropped_gses <- setdiff(all_gses, comb_gses)
    dropped_str <- if (length(dropped_gses) == 0) "None" else paste(dropped_gses, collapse = "; ")

    if (!is.null(comb_dir) && !dir.exists(comb_dir)) {
      dir.create(comb_dir, recursive = TRUE)
    }

    comb_comps <- names(comparison_gse_map)[comparison_gse_map %in% comb_gses]
    sub_lfc <- lfc_matrix[, comb_comps, drop = FALSE]
    sub_se <- se_matrix[, comb_comps, drop = FALSE]
    sub_pval <- pval_matrix[, comb_comps, drop = FALSE]
    sub_gse_map <- comparison_gse_map[comb_comps]

    # Check for existing completed results to reuse
    results_file <- if (!is.null(comb_dir)) file.path(comb_dir, "DESeq2_contrast_results.csv") else NULL
    sig_3lre_file <- if (!is.null(comb_dir)) file.path(comb_dir, "DESeq2_contrast_results_significant_Meta3LRE.csv") else NULL
    sig_fe_file <- if (!is.null(comb_dir)) file.path(comb_dir, "DESeq2_contrast_results_significant_MetaFE.csv") else NULL
    info_file <- if (!is.null(comb_dir)) file.path(comb_dir, "combination_info.txt") else NULL

    files_complete <- !is.null(comb_dir) && file.exists(results_file) && file.exists(sig_3lre_file) && file.exists(sig_fe_file) && file.exists(info_file)

    if (!force_rerun && files_complete) {
      info_lines <- tryCatch(readLines(info_file), error = function(e) character())
      inc_line <- grep("^Included GSE Datasets:", info_lines, value = TRUE)
      m_line <- grep("^Min Datasets Flag \\(-m\\):", info_lines, value = TRUE)
      expected_inc_str <- paste0("Included GSE Datasets: ", paste(comb_gses, collapse = ", "))
      expected_m_str <- paste0("Min Datasets Flag (-m): ", min_datasets)

      datasets_match <- length(inc_line) > 0 && trimws(inc_line[1]) == trimws(expected_inc_str)
      m_param_match <- length(m_line) == 0 || trimws(m_line[1]) == trimws(expected_m_str)

      if (datasets_match && m_param_match) {
        res_cached <- tryCatch(read.csv(results_file, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(res_cached) && "Meta_3LRE_padj" %in% colnames(res_cached) && "Meta_FE_padj" %in% colnames(res_cached)) {
          cat(sprintf("  ✓ Found completed results for %s. Reusing cached results without re-fitting.\n", comb_name))

          sig_3lre <- !is.na(res_cached$Meta_3LRE_padj) & res_cached$Meta_3LRE_padj < 0.05 & abs(res_cached$Meta_3LRE_LFC) > 0.5
          n_3lre_tot <- sum(sig_3lre)
          n_3lre_up <- sum(sig_3lre & res_cached$Meta_3LRE_LFC > 0.5)
          n_3lre_down <- sum(sig_3lre & res_cached$Meta_3LRE_LFC < -0.5)

          sig_fe <- !is.na(res_cached$Meta_FE_padj) & res_cached$Meta_FE_padj < 0.05 & abs(res_cached$Meta_FE_LFC) > 0.5
          n_fe_tot <- sum(sig_fe)
          n_fe_up <- sum(sig_fe & res_cached$Meta_FE_LFC > 0.5)
          n_fe_down <- sum(sig_fe & res_cached$Meta_FE_LFC < -0.5)

          n_analyzed <- sum(!is.na(res_cached$Meta_3LRE_pval) | !is.na(res_cached$Meta_FE_pval))

          summary_row <- data.frame(
            Combination = comb_name,
            K_datasets = k_cur,
            Included_GSEs = paste(comb_gses, collapse = "; "),
            Dropped_GSEs = dropped_str,
            Num_Comparisons = length(comb_comps),
            Num_Analyzed_Genes = n_analyzed,
            Meta3LRE_Total_DEG = n_3lre_tot,
            Meta3LRE_Up_DEG = n_3lre_up,
            Meta3LRE_Down_DEG = n_3lre_down,
            MetaFE_Total_DEG = n_fe_tot,
            MetaFE_Up_DEG = n_fe_up,
            MetaFE_Down_DEG = n_fe_down,
            stringsAsFactors = FALSE
          )

          cat(sprintf("  ✓ Loaded %s: Analyzed %d genes | 3LRE DEGs: %d (Up: %d, Down: %d) | FE DEGs: %d (Up: %d, Down: %d)\n\n",
                      comb_name, n_analyzed, n_3lre_tot, n_3lre_up, n_3lre_down, n_fe_tot, n_fe_up, n_fe_down))
          return(summary_row)
        }
      }
    }

    # Effective min datasets for gene inclusion in this combination
    eff_min_ds <- if (min_datasets == 0) k_cur else min(k_cur, min_datasets)
    eff_min_ds <- max(2, eff_min_ds)

    valid_mat <- !is.na(sub_lfc) & !is.na(sub_se) & (sub_se > 0) & !is.na(sub_pval)

    # Determine dataset count per gene
    gse_has_gene <- matrix(FALSE, nrow = n_genes, ncol = length(comb_gses))
    colnames(gse_has_gene) <- comb_gses
    for (gse in comb_gses) {
      comps_gse <- comb_comps[sub_gse_map == gse]
      if (length(comps_gse) == 1) {
        gse_has_gene[, gse] <- valid_mat[, comps_gse]
      } else if (length(comps_gse) > 1) {
        gse_has_gene[, gse] <- rowSums(valid_mat[, comps_gse, drop = FALSE]) > 0
      }
    }
    ds_count_per_gene <- rowSums(gse_has_gene)
    genes_to_fit <- which(ds_count_per_gene >= eff_min_ds)

    cat(sprintf("  - %d / %d genes qualify (>= %d datasets). Running meta-analysis...\n",
                length(genes_to_fit), n_genes, eff_min_ds))

    # Output vectors
    re3L_lfc_vec <- rep(NA_real_, n_genes)
    re3L_se_vec <- rep(NA_real_, n_genes)
    re3L_pval_vec <- rep(NA_real_, n_genes)
    fe_lfc_vec <- rep(NA_real_, n_genes)
    fe_se_vec <- rep(NA_real_, n_genes)
    fe_pval_vec <- rep(NA_real_, n_genes)
    tau2_L3_vec <- rep(NA_real_, n_genes)
    tau2_L2_vec <- rep(NA_real_, n_genes)
    I2_L3_vec <- rep(NA_real_, n_genes)
    I2_L2_vec <- rep(NA_real_, n_genes)
    het_pval_vec <- rep(NA_real_, n_genes)

    if (length(genes_to_fit) > 0) {
      pb <- txtProgressBar(min = 0, max = length(genes_to_fit), style = 3)
      for (step_i in seq_along(genes_to_fit)) {
        setTxtProgressBar(pb, step_i)
        g_idx <- genes_to_fit[step_i]
        gene_id <- all_genes[g_idx]

        val_cols <- which(valid_mat[g_idx, ])
        yi_val <- sub_lfc[g_idx, val_cols]
        sei_val <- sub_se[g_idx, val_cols]
        vi_val <- sei_val^2
        val_comps <- comb_comps[val_cols]
        val_gse <- sub_gse_map[val_cols]
        n_val <- length(val_comps)

        # Construct variance-covariance matrix
        V_exact <- diag(vi_val, nrow = n_val, ncol = n_val)
        rownames(V_exact) <- val_comps
        colnames(V_exact) <- val_comps

        for (r in 1:n_val) {
          comp1 <- val_comps[r]
          gse1 <- val_gse[r]
          for (c in 1:n_val) {
            comp2 <- val_comps[c]
            gse2 <- val_gse[c]
            if (r != c && gse1 == gse2) {
              true_cov <- NA
              if (!is.null(all_cov_mats[[gse1]])) {
                if ("dataset_level" %in% names(all_cov_mats[[gse1]])) {
                  true_cov <- all_cov_mats[[gse1]][["dataset_level"]][comp1, comp2]
                } else if (!is.null(all_cov_mats[[gse1]][[gene_id]])) {
                  true_cov <- all_cov_mats[[gse1]][[gene_id]][comp1, comp2]
                }
              }
              if (!is.na(true_cov)) {
                s1 <- sign_map[[comp1]]; if (is.null(s1) || is.na(s1)) s1 <- 1
                s2 <- sign_map[[comp2]]; if (is.null(s2) || is.na(s2)) s2 <- 1
                V_exact[r, c] <- true_cov * s1 * s2
              }
            }
          }
        }

        V_exact <- (V_exact + t(V_exact)) / 2
        is_pd <- tryCatch({ chol(V_exact); TRUE }, error = function(e) FALSE)
        if (!is_pd) {
          diag(V_exact) <- diag(V_exact) + 1e-6
        }

        has_mult <- any(table(val_gse) > 1)
        rnd_f <- if (has_mult) ~ 1 | val_gse / val_comps else ~ 1 | val_gse

        # 3LRE Model
        res_3lre <- tryCatch(
          metafor::rma.mv(yi = yi_val, V = V_exact, random = rnd_f, method = "REML"),
          error = function(e) {
            tryCatch(
              metafor::rma.mv(yi = yi_val, V = V_exact, random = rnd_f, method = "REML",
                              control = list(optimizer = "optim", optmethod = "Nelder-Mead")),
              error = function(e2) NULL
            )
          }
        )

        # FE Model
        res_fe <- tryCatch(
          metafor::rma.mv(yi = yi_val, V = V_exact),
          error = function(e) NULL
        )

        # 3LRE Extraction
        if (!is.null(res_3lre)) {
          rob_res <- tryCatch(
            metafor::robust(res_3lre, cluster = val_gse, clubSandwich = TRUE),
            error = function(e) NULL
          )
          if (!is.null(rob_res)) {
            re3L_lfc_vec[g_idx] <- as.numeric(rob_res$beta)[1]
            re3L_se_vec[g_idx] <- as.numeric(rob_res$se)[1]
            re3L_pval_vec[g_idx] <- as.numeric(rob_res$pval)[1]
          }
          if (length(res_3lre$sigma2) >= 2) {
            tau2_L3_vec[g_idx] <- res_3lre$sigma2[1]
            tau2_L2_vec[g_idx] <- res_3lre$sigma2[2]
            i2_res <- calculate_pseudo_i2_dropout(vi_val, c(res_3lre$sigma2[1], res_3lre$sigma2[2]))
            I2_L3_vec[g_idx] <- i2_res[1]
            I2_L2_vec[g_idx] <- i2_res[2]
          } else {
            tau2_L3_vec[g_idx] <- res_3lre$sigma2[1]
            tau2_L2_vec[g_idx] <- 0
            i2_res <- calculate_pseudo_i2_dropout(vi_val, c(res_3lre$sigma2[1]))
            I2_L3_vec[g_idx] <- i2_res[1]
            I2_L2_vec[g_idx] <- 0
          }
          het_pval_vec[g_idx] <- res_3lre$QEp
        }

        # FE Extraction
        if (!is.null(res_fe)) {
          fe_lfc_vec[g_idx] <- as.numeric(res_fe$beta)[1]
          fe_se_vec[g_idx] <- as.numeric(res_fe$se)[1]
          fe_pval_vec[g_idx] <- res_fe$pval[1]
        }
      }
      close(pb)
    }

    # BH FDR adjustments
    padj_3lre <- p.adjust(re3L_pval_vec, method = "BH")
    padj_fe <- p.adjust(fe_pval_vec, method = "BH")
    padj_het <- p.adjust(het_pval_vec, method = "BH")

    # Assemble result data frame
    res_comb_df <- data.frame(
      cluster_id = all_genes,
      Meta_3LRE_LFC = re3L_lfc_vec,
      Meta_3LRE_SE = re3L_se_vec,
      Meta_3LRE_pval = re3L_pval_vec,
      Meta_3LRE_padj = padj_3lre,
      Meta_FE_LFC = fe_lfc_vec,
      Meta_FE_SE = fe_se_vec,
      Meta_FE_pval = fe_pval_vec,
      Meta_FE_padj = padj_fe,
      Heterogeneity_3L_I2_Level3 = I2_L3_vec,
      Heterogeneity_3L_I2_Level2 = I2_L2_vec,
      Heterogeneity_3L_tau2_Level3 = tau2_L3_vec,
      Heterogeneity_3L_tau2_Level2 = tau2_L2_vec,
      Heterogeneity_3L_pval = het_pval_vec,
      Heterogeneity_3L_padj = padj_het,
      stringsAsFactors = FALSE
    )

    res_comb_export <- left_join(data.frame(cluster_id = all_genes, stringsAsFactors = FALSE), cluster_annotation, by = "cluster_id")
    res_comb_export <- left_join(res_comb_export, res_comb_df, by = "cluster_id")

    # Significant subsets
    res_comb_sig_3lre <- res_comb_export %>%
      dplyr::filter(!is.na(Meta_3LRE_padj) & Meta_3LRE_padj < 0.05 & abs(Meta_3LRE_LFC) > 0.5)

    res_comb_sig_fe <- res_comb_export %>%
      dplyr::filter(!is.na(Meta_FE_padj) & Meta_FE_padj < 0.05 & abs(Meta_FE_LFC) > 0.5)

    # Counts
    n_3lre_tot <- nrow(res_comb_sig_3lre)
    n_3lre_up <- sum(!is.na(res_comb_sig_3lre$Meta_3LRE_LFC) & res_comb_sig_3lre$Meta_3LRE_LFC > 0.5)
    n_3lre_down <- sum(!is.na(res_comb_sig_3lre$Meta_3LRE_LFC) & res_comb_sig_3lre$Meta_3LRE_LFC < -0.5)

    n_fe_tot <- nrow(res_comb_sig_fe)
    n_fe_up <- sum(!is.na(res_comb_sig_fe$Meta_FE_LFC) & res_comb_sig_fe$Meta_FE_LFC > 0.5)
    n_fe_down <- sum(!is.na(res_comb_sig_fe$Meta_FE_LFC) & res_comb_sig_fe$Meta_FE_LFC < -0.5)

    n_analyzed <- sum(!is.na(re3L_pval_vec) | !is.na(fe_pval_vec))

    # Export CSVs and metadata for this combination if comb_dir provided
    if (!is.null(comb_dir)) {
      write.csv(res_comb_export, file.path(comb_dir, "DESeq2_contrast_results.csv"), row.names = FALSE)
      write.csv(res_comb_sig_3lre, file.path(comb_dir, "DESeq2_contrast_results_significant_Meta3LRE.csv"), row.names = FALSE)
      write.csv(res_comb_sig_fe, file.path(comb_dir, "DESeq2_contrast_results_significant_MetaFE.csv"), row.names = FALSE)

      # Write detailed combination metadata file
      info_lines <- c(
        paste0("Combination: ", comb_name),
        paste0("Total GSEs in Project: ", N),
        paste0("Number of Included GSEs (K): ", k_cur),
        paste0("Included GSE Datasets: ", paste(comb_gses, collapse = ", ")),
        paste0("Excluded/Dropped GSE Datasets: ", dropped_str),
        paste0("Min Datasets Flag (-m): ", min_datasets),
        paste0("Number of Included Comparisons: ", length(comb_comps)),
        paste0("Included Comparisons: ", paste(comb_comps, collapse = ", ")),
        paste0("Number of Analyzed Genes: ", n_analyzed),
        paste0("Meta3LRE DEGs (FDR < 0.05, |LFC| > 0.5): ", n_3lre_tot, " (Up: ", n_3lre_up, ", Down: ", n_3lre_down, ")"),
        paste0("MetaFE DEGs (FDR < 0.05, |LFC| > 0.5): ", n_fe_tot, " (Up: ", n_fe_up, ", Down: ", n_fe_down, ")")
      )
      writeLines(info_lines, file.path(comb_dir, "combination_info.txt"))
    }

    summary_row <- data.frame(
      Combination = comb_name,
      K_datasets = k_cur,
      Included_GSEs = paste(comb_gses, collapse = "; "),
      Dropped_GSEs = dropped_str,
      Num_Comparisons = length(comb_comps),
      Num_Analyzed_Genes = n_analyzed,
      Meta3LRE_Total_DEG = n_3lre_tot,
      Meta3LRE_Up_DEG = n_3lre_up,
      Meta3LRE_Down_DEG = n_3lre_down,
      MetaFE_Total_DEG = n_fe_tot,
      MetaFE_Up_DEG = n_fe_up,
      MetaFE_Down_DEG = n_fe_down,
      stringsAsFactors = FALSE
    )

    cat(sprintf("  ✓ Finished %s: Analyzed %d genes | 3LRE DEGs: %d (Up: %d, Down: %d) | FE DEGs: %d (Up: %d, Down: %d)\n\n",
                comb_name, n_analyzed, n_3lre_tot, n_3lre_up, n_3lre_down, n_fe_tot, n_fe_up, n_fe_down))

    return(summary_row)
  }

  summary_rows <- list()

  # 1. Baseline Row (Full N-Dataset Reference)
  cat("[Baseline] Evaluating Full Baseline (All datasets) reference (3LRE & FE)...\n")
  base_dir <- file.path(dropout_dir, "baseline")
  summary_rows[[1]] <- fit_one_combination(all_gses, "Baseline_All", base_dir)

  # 2. Iterate over drop-out combinations
  for (i in seq_along(comb_list)) {
    comb_gses <- comb_list[[i]]
    comb_name <- paste0("combination_", i)
    comb_dir <- file.path(dropout_dir, comb_name)

    cat(sprintf("[%d/%d] %s (K=%d): Included: [%s] | Dropped: [%s]\n",
                i, n_combs, comb_name, length(comb_gses),
                paste(comb_gses, collapse = ", "),
                paste(setdiff(all_gses, comb_gses), collapse = ", ")))

    summary_rows[[length(summary_rows) + 1]] <- fit_one_combination(comb_gses, comb_name, comb_dir)
  }

  summary_df <- bind_rows(summary_rows)
  summary_file <- file.path(dropout_dir, "dropout_summary.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)

  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("✓ DATASET DROP-OUT TEST COMPLETE!\n")
  cat("Master summary saved to:", summary_file, "\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  print(summary_df[, c("Combination", "K_datasets", "Dropped_GSEs", "Num_Analyzed_Genes", "Meta3LRE_Total_DEG", "MetaFE_Total_DEG")])
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  return(invisible(summary_df))
}
