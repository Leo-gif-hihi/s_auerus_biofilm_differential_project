# ==============================================================================
# Script: imodulon_decoupler_analysis.R
# Purpose: Infers upstream transcription factor and iModulon regulatory activity
#          using decoupleR via two complementary approaches:
#            Approach 1: Directly on pooled meta-analysis gene statistics
#                        (ULM + Consensus across Meta_3LRE, Meta_FE, DirectionalStouffer)
#            Approach 2: Per-dataset ULM inference followed by 3-level
#                        random-effects & fixed-effects meta-analysis on Fisher-Z
#                        (mirroring ora_meta_analysis.R)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(decoupleR)
  library(metafor)
  library(clubSandwich)
})

if (!exists("map_clusters_to_imodulon") && file.exists("gene_id_mapping.R")) {
  source("gene_id_mapping.R")
}

#' Curate STAPH_PRECISE165 iModulons and build weighted network
#' @keywords internal
.curate_imodulon_network <- function(m_file, table_file, thresholds_file, weight_threshold = 0.08) {
  cat("Loading iModulonDB reference files...\n")
  M_mat <- read.csv(m_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  imod_table <- read.csv(table_file, stringsAsFactors = FALSE)

  # Exclude phage/mobile-element artifacts and single-gene ICA components:
  # - Prophages (4): Phi-Sa3-1/2/3, Genomic Island 3
  # - Single Gene (8): SG_1..SG_8
  # - Technical artifacts: scfm2 artefact
  # - Pure genomic island modules: SCCMec and ACME (category="Other", regulator="Genomic Island 1")
  # Keeps Uncharacterized-1..9 and SaeR-2 (category="Virulence", regulator="SaeR, Genomic Island 8")
  excluded_cats <- c("Prophages", "Genomic Islands", "Single Gene")
  clean_table <- imod_table %>%
    filter(
      !category %in% excluded_cats,
      !grepl("artefact|artifact", name, ignore.case = TRUE),
      !grepl("^SG_", name),
      !grepl("^Genomic Island", ifelse(is.na(regulator), "", regulator), ignore.case = TRUE)
    )

  dropped <- setdiff(colnames(M_mat), clean_table$name)
  valid_imodulons <- intersect(colnames(M_mat), clean_table$name)
  cat(sprintf(
    "  Retained %d / %d iModulons (excluded %d: %s).\n",
    length(valid_imodulons), ncol(M_mat), length(dropped), paste(dropped, collapse = ", ")
  ))

  M_mat <- M_mat[, valid_imodulons, drop = FALSE]

  # Resolve per-iModulon weight cutoffs (official thresholds preferred, scalar fallback)
  use_official <- !is.null(thresholds_file) && file.exists(thresholds_file)
  if (use_official) {
    thr_df <- read.csv(thresholds_file, stringsAsFactors = FALSE)
    thr_map <- setNames(as.numeric(thr_df$threshold), thr_df$imodulon_name)
    missing_thr <- setdiff(valid_imodulons, names(thr_map))
    if (length(missing_thr) > 0) {
      thr_map[missing_thr] <- weight_threshold
    }
    cat(sprintf("  Using official per-iModulon thresholds from %s.\n", thresholds_file))
  } else {
    thr_map <- setNames(rep(weight_threshold, length(valid_imodulons)), valid_imodulons)
    cat(sprintf("  Thresholds file not found; using fixed |weight| >= %.3f.\n", weight_threshold))
  }

  net_raw <- M_mat %>%
    rownames_to_column(var = "target") %>%
    pivot_longer(cols = -target, names_to = "source", values_to = "weight") %>%
    mutate(cutoff = unname(thr_map[source])) %>%
    filter(abs(weight) >= cutoff) %>%
    mutate(mor = weight) %>%
    select(source, target, mor)

  clean_table <- clean_table %>%
    mutate(
      has_reg = !is.na(regulator) & trimws(regulator) != "" & trimws(regulator) != "-",
      source_annot = ifelse(has_reg, paste0(name, " (", regulator, ")"), name)
    )

  annot_map <- setNames(clean_table$source_annot, clean_table$name)
  net_raw$source_annot <- unname(annot_map[net_raw$source])

  imod_meta <- clean_table %>%
    filter(name %in% valid_imodulons) %>%
    transmute(
      source = source_annot,
      imodulon_name = name,
      category = category,
      regulator = ifelse(has_reg, regulator, NA_character_),
      explained_variance = explained_variance,
      n_genes_imodulon = num_genes
    )

  list(
    m_gene_ids = rownames(M_mat),
    net_raw = net_raw,
    imod_meta = imod_meta
  )
}

#' Build cluster_id-mapped decoupleR network and export mapping diagnostics
#' @keywords internal
.build_mapped_network <- function(net_raw,
                                  cluster_df,
                                  m_gene_ids,
                                  gene_info_file,
                                  trn_file,
                                  orthologue_file,
                                  annote_dir,
                                  operon_file,
                                  gaf_file,
                                  out_dir) {
  cat("Harmonizing pipeline cluster_ids with iModulon target genes...\n")
  map_res <- map_clusters_to_imodulon(
    cluster_df = cluster_df,
    m_gene_ids = m_gene_ids,
    gene_info_file = gene_info_file,
    trn_file = trn_file,
    orthologue_file = orthologue_file,
    annote_dir = annote_dir,
    operon_file = operon_file,
    gaf_file = gaf_file,
    verbose = TRUE
  )

  if (nrow(map_res$mapping_df) == 0) {
    stop("CRITICAL ERROR: Zero overlapping genes between pipeline cluster IDs and iModulon M matrix.")
  }

  net_mapped <- net_raw %>%
    inner_join(map_res$mapping_df, by = "target", relationship = "many-to-many") %>%
    group_by(source = source_annot, target = cluster_id) %>%
    summarise(
      mor = mean(mor),
      n_genes_collapsed = n(),
      .groups = "drop"
    )

  if (nrow(net_mapped) == 0) {
    stop("CRITICAL ERROR: Zero overlapping interactions between pipeline cluster IDs and thresholded iModulons.")
  }

  n_multi <- sum(net_mapped$n_genes_collapsed > 1)
  if (n_multi > 0) {
    cat(sprintf("  ℹ Aggregated %d cluster-iModulon pairs spanning multiple mapped target tags (mean mor).\n", n_multi))
  }

  n_net_clusters <- length(unique(net_mapped$target))
  n_net_targets <- length(unique(net_raw$target))
  n_raw_targets_mapped <- length(unique(intersect(net_raw$target, map_res$mapping_df$target)))
  pct_net_clusters <- 100 * n_net_clusters / map_res$stats$n_total_clusters
  pct_net_targets <- if (n_net_targets > 0) 100 * n_raw_targets_mapped / n_net_targets else 0

  cat(sprintf(
    "  ✓ Thresholded network: %d edges | %d / %d cluster_ids in network (%.2f%%) | %d / %d thresholded M-genes mapped (%.2f%%).\n",
    nrow(net_mapped), n_net_clusters, map_res$stats$n_total_clusters, pct_net_clusters,
    n_raw_targets_mapped, n_net_targets, pct_net_targets
  ))

  # Save mapping diagnostics
  summary_tbl <- data.frame(
    Metric = c(
      "Total_Pipeline_Clusters",
      "Clusters_Mapped_to_M_Matrix",
      "Pct_Clusters_Mapped_to_M_Matrix",
      "Clusters_Mapped_to_Multiple_Targets",
      "Pct_Clusters_Mapped_to_Multiple_Targets",
      "Clusters_in_Thresholded_Network",
      "Pct_Clusters_in_Thresholded_Network",
      "Total_M_Matrix_Genes",
      "M_Matrix_Genes_Mapped",
      "Pct_M_Matrix_Genes_Mapped",
      "Thresholded_Network_Genes",
      "Thresholded_Network_Genes_Mapped",
      "Pct_Thresholded_Network_Genes_Mapped",
      "Total_Network_Edges"
    ),
    Value = c(
      map_res$stats$n_total_clusters,
      map_res$stats$n_mapped_clusters,
      round(map_res$stats$pct_mapped_clusters, 2),
      map_res$stats$n_multi_target_clusters,
      round(map_res$stats$pct_multi_target_clusters, 2),
      n_net_clusters,
      round(pct_net_clusters, 2),
      map_res$stats$n_total_m_genes,
      map_res$stats$n_mapped_m_genes,
      round(map_res$stats$pct_mapped_m_genes, 2),
      n_net_targets,
      n_raw_targets_mapped,
      round(pct_net_targets, 2),
      nrow(net_mapped)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(summary_tbl, file.path(out_dir, "iModulon_Gene_Mapping_Summary.csv"), row.names = FALSE)
  write.csv(map_res$cluster_summary, file.path(out_dir, "iModulon_Cluster_Mapping.csv"), row.names = FALSE)

  list(net_mapped = net_mapped, mapping_stats = map_res$stats)
}

#' Run decoupleR safely (handles MLM collinearity fallback if needed)
#' @keywords internal
.run_decoupler_safe <- function(input_mat, net_mapped, min_targets = 5, run_consensus = TRUE, wmean_times = 500) {
  if (!run_consensus) {
    return(decoupleR::run_ulm(
      mat = input_mat,
      network = net_mapped,
      .source = source,
      .target = target,
      .mor = mor,
      minsize = min_targets
    ))
  }

  tryCatch(
    decoupleR::decouple(
      mat = input_mat,
      network = net_mapped,
      .source = source,
      .target = target,
      statistics = c("ulm", "wmean", "mlm"),
      args = list(wmean = list(times = wmean_times, seed = 42)),
      consensus_score = TRUE,
      minsize = min_targets
    ),
    error = function(e) {
      warning(sprintf("decouple() with MLM failed (%s); falling back to ULM + WMEAN.", e$message))
      decoupleR::decouple(
        mat = input_mat,
        network = net_mapped,
        .source = source,
        .target = target,
        statistics = c("ulm", "wmean"),
        args = list(wmean = list(times = wmean_times, seed = 42)),
        consensus_score = TRUE,
        minsize = min_targets
      )
    }
  )
}

#' Generate horizontal barplot of top significant iModulon drivers
#' @keywords internal
.plot_imodulon_drivers <- function(df, score_col, p_col, padj_col, title_str, subtitle_str, x_lab, out_file) {
  sub_sig <- df %>% filter(!is.na(.data[[padj_col]]) & .data[[padj_col]] < 0.05)
  sig_note <- "BH padj < 0.05"
  if (nrow(sub_sig) == 0) {
    sub_sig <- df %>% filter(!is.na(.data[[p_col]]) & .data[[p_col]] < 0.05)
    sig_note <- "nominal p < 0.05 (none BH padj < 0.05)"
  }

  top_hits <- sub_sig %>%
    arrange(desc(abs(.data[[score_col]]))) %>%
    head(20)

  if (nrow(top_hits) == 0) {
    cat(sprintf("  ⚠ No iModulons met p < 0.05 for plot: %s\n", basename(out_file)))
    return(invisible(NULL))
  }

  top_hits$is_pos <- top_hits[[score_col]] > 0
  p <- ggplot(top_hits, aes(x = reorder(source, .data[[score_col]]), y = .data[[score_col]], fill = is_pos)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(
      values = c("TRUE" = "#D95F02", "FALSE" = "#1F78B4"),
      labels = c("TRUE" = "Activated in Biofilm", "FALSE" = "Repressed in Biofilm"),
      drop = FALSE
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10),
      legend.position = "bottom"
    ) +
    labs(
      title = title_str,
      subtitle = paste0(subtitle_str, " | ", sig_note),
      x = "iModulon / Regulator",
      y = x_lab,
      fill = "State"
    )

  ggsave(out_file, plot = p, width = 9.5, height = 7, dpi = 300)
  cat("  ✓ Saved plot:", out_file, "\n")
}

#' Compute effective network target counts present in a gene vector
#' @keywords internal
.count_mapped_targets_per_source <- function(net_mapped, measured_ids, min_targets) {
  net_sub <- net_mapped[net_mapped$target %in% measured_ids, , drop = FALSE]
  cnts <- table(net_sub$source)
  cnts <- cnts[cnts >= min_targets]
  setNames(as.integer(cnts), names(cnts))
}

#' Fit 3-level random-effects model with REML + Nelder-Mead fallback + clubSandwich RVE
#' @keywords internal
.fit_3lre_robust <- function(yi, V_mat, gse_vec, comp_vec, W_vec = NULL) {
  fit_fn <- function(ctrl = list()) {
    if (is.null(W_vec)) {
      metafor::rma.mv(yi = yi, V = V_mat, random = ~ 1 | gse_vec / comp_vec, method = "REML", control = ctrl)
    } else {
      metafor::rma.mv(yi = yi, V = V_mat, W = W_vec, random = ~ 1 | gse_vec / comp_vec, method = "REML", control = ctrl)
    }
  }
  res_re <- tryCatch(fit_fn(), error = function(e) {
    tryCatch(fit_fn(list(optimizer = "optim", optmethod = "Nelder-Mead")), error = function(e2) NULL)
  })
  if (is.null(res_re)) return(NULL)

  rob <- tryCatch(metafor::robust(res_re, cluster = gse_vec, clubSandwich = TRUE), error = function(e) NULL)
  list(re = res_re, rob = rob)
}

#' Calculate pseudo-I2 for 3-level meta-analysis (mirrors DESeq2 main script)
#' @keywords internal
.calc_pseudo_i2 <- function(vi, sigma2) {
  w <- 1 / vi
  s <- sum(w) - sum(w^2) / sum(w)
  if (!is.finite(s) || s <= 0) return(c(NA_real_, NA_real_))
  v_typ <- (length(vi) - 1) / s
  tot <- sum(sigma2) + v_typ
  if (!is.finite(tot) || tot <= 0) return(c(NA_real_, NA_real_))
  100 * sigma2 / tot
}

#' Run Approach 1: decoupleR directly on meta-analyzed gene statistics
#' @keywords internal
.run_approach1_direct <- function(res_df, net_mapped, imod_meta, sample_label, out_dir, vis_dir, min_targets, run_consensus) {
  cat("\n--- Approach 1: Direct decoupleR on Pooled Meta-Analysis Gene Statistics ---\n")

  stat_df <- res_df %>%
    mutate(
      Meta_3LRE_Stat = ifelse(!is.na(Meta_3LRE_SE) & Meta_3LRE_SE > 0 & !is.na(Meta_3LRE_LFC),
        Meta_3LRE_LFC / Meta_3LRE_SE, NA_real_
      ),
      Meta_FE_Stat = ifelse(!is.na(Meta_FE_SE) & Meta_FE_SE > 0 & !is.na(Meta_FE_LFC),
        Meta_FE_LFC / Meta_FE_SE, NA_real_
      ),
      Meta_DirectionalStouffer_Stat = if ("Meta_DirectionalStouffer_Z" %in% colnames(.)) {
        ifelse(!is.na(Meta_DirectionalStouffer_Z), Meta_DirectionalStouffer_Z, NA_real_)
      } else {
        NA_real_
      }
    )

  model_specs <- list(
    Meta3LRE = "Meta_3LRE_Stat",
    MetaFE = "Meta_FE_Stat",
    MetaDirectionalStouffer = "Meta_DirectionalStouffer_Stat"
  )

  direct_results <- list()

  for (model_label in names(model_specs)) {
    stat_col <- model_specs[[model_label]]
    valid_df <- stat_df[!is.na(stat_df[[stat_col]]) & is.finite(stat_df[[stat_col]]), , drop = FALSE]
    if (nrow(valid_df) < min_targets) {
      cat(sprintf("  ⚠ Skipping [%s]: fewer than %d valid genes.\n", model_label, min_targets))
      next
    }

    input_mat <- matrix(valid_df[[stat_col]], ncol = 1, dimnames = list(valid_df$cluster_id, sample_label))
    cat(sprintf("  Running decoupleR [%s] (col=%s, N=%d genes)...\n", model_label, sample_label, nrow(input_mat)))

    scores <- .run_decoupler_safe(input_mat, net_mapped, min_targets = min_targets, run_consensus = run_consensus, wmean_times = 1000)
    target_counts <- .count_mapped_targets_per_source(net_mapped, rownames(input_mat), min_targets)

    export_df <- scores %>%
      group_by(statistic) %>%
      mutate(padj = p.adjust(p_value, method = "BH")) %>%
      ungroup() %>%
      mutate(n_targets_mapped = unname(target_counts[source])) %>%
      left_join(imod_meta, by = "source") %>%
      arrange(statistic, desc(abs(score)))

    csv_file <- file.path(out_dir, paste0("iModulon_activities_", model_label, ".csv"))
    write.csv(export_df, csv_file, row.names = FALSE)
    cat("  ✓ Exported:", csv_file, "\n")

    if (run_consensus) {
      .plot_imodulon_drivers(
        df = export_df %>% filter(statistic == "consensus"),
        score_col = "score",
        p_col = "p_value",
        padj_col = "padj",
        title_str = paste("Top iModulon Drivers (Direct):", model_label, "[", sample_label, "]"),
        subtitle_str = "Method: Consensus (ULM/WMEAN/MLM)",
        x_lab = "Inferred Activity Score (Consensus Z)",
        out_file = file.path(vis_dir, paste0("Barplot_iModulon_Drivers_", model_label, ".png"))
      )
    }

    .plot_imodulon_drivers(
      df = export_df %>% filter(statistic == "ulm"),
      score_col = "score",
      p_col = "p_value",
      padj_col = "padj",
      title_str = paste("Top iModulon Drivers (Direct ULM):", model_label, "[", sample_label, "]"),
      subtitle_str = "Method: ULM (Univariate Linear Model)",
      x_lab = "Inferred Activity Score (ULM t-stat)",
      out_file = file.path(vis_dir, paste0("Barplot_iModulon_Drivers_", model_label, "_ULM.png"))
    )

    direct_results[[model_label]] <- export_df
  }

  direct_results
}

#' Run Approach 2: Per-dataset ULM -> 3LRE + FE Meta-Analysis (mirrors ora_meta_analysis.R)
#' @keywords internal
.run_approach2_per_dataset_meta <- function(lfc_list,
                                            se_list,
                                            comparison_gse_map,
                                            comparison_control_map,
                                            neff_list = NULL,
                                            min_ds,
                                            net_mapped,
                                            imod_meta,
                                            sample_label,
                                            out_dir,
                                            vis_dir,
                                            min_targets) {
  cat("\n--- Approach 2: Per-Dataset ULM -> 3LRE / FE Meta-Analysis ---\n")
  if (is.null(lfc_list) || is.null(se_list) || length(lfc_list) == 0) {
    cat("  ⚠ lfc_list / se_list not provided. Skipping Approach 2.\n")
    return(NULL)
  }

  comp_names <- intersect(names(lfc_list), names(se_list))
  if (length(comp_names) == 0) {
    cat("  ⚠ No overlapping comparisons between lfc_list and se_list. Skipping Approach 2.\n")
    return(NULL)
  }

  per_comp_scores <- list()
  comp_n_measured <- setNames(integer(length(comp_names)), comp_names)
  comp_target_counts <- list()

  for (comp_name in comp_names) {
    lfcs <- lfc_list[[comp_name]]
    ses <- se_list[[comp_name]]
    if (!is.null(names(lfcs)) && !is.null(names(ses))) {
      common_ids <- intersect(names(lfcs), names(ses))
      lfcs <- lfcs[common_ids]
      ses <- ses[common_ids]
    }
    valid_mask <- !is.na(lfcs) & !is.na(ses) & is.finite(lfcs) & is.finite(ses) & ses > 0
    if (sum(valid_mask) < max(min_targets, 10L)) next

    stat_vec <- lfcs[valid_mask] / ses[valid_mask]
    comp_mat <- matrix(stat_vec, ncol = 1, dimnames = list(names(stat_vec), comp_name))
    comp_n_measured[comp_name] <- nrow(comp_mat)
    comp_target_counts[[comp_name]] <- .count_mapped_targets_per_source(net_mapped, rownames(comp_mat), min_targets)

    sc <- .run_decoupler_safe(comp_mat, net_mapped, min_targets = min_targets, run_consensus = FALSE)
    sc$gse <- if (!is.null(comparison_gse_map) && comp_name %in% names(comparison_gse_map)) {
      unname(comparison_gse_map[comp_name])
    } else {
      sub("_.*", "", comp_name)
    }
    sc$control <- if (!is.null(comparison_control_map) && comp_name %in% names(comparison_control_map)) {
      unname(comparison_control_map[comp_name])
    } else {
      NA_character_
    }
    sc$N_measured <- nrow(comp_mat)
    per_comp_scores[[comp_name]] <- sc
  }

  if (length(per_comp_scores) == 0) {
    cat("  ⚠ No comparisons produced valid decoupleR scores.\n")
    return(NULL)
  }

  all_per_df <- bind_rows(per_comp_scores) %>%
    group_by(condition, statistic) %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    left_join(imod_meta %>% select(source, imodulon_name, category, regulator), by = "source")

  per_csv <- file.path(out_dir, "iModulon_PerDataset_Activities.csv")
  write.csv(all_per_df, per_csv, row.names = FALSE)
  cat(sprintf("  ✓ Exported per-dataset ULM activities (%d comparisons) -> %s\n", length(per_comp_scores), per_csv))

  # --- Pool each iModulon (ULM Fisher-Z) across comparisons ---
  sources <- unique(all_per_df$source)
  min_ds_req <- max(2L, as.integer(min_ds))
  ulm_df <- all_per_df %>% filter(statistic == "ulm")

  meta_rows <- vector("list", length(sources))
  names(meta_rows) <- sources

  for (i in seq_along(sources)) {
    src <- sources[i]
    u_sub <- ulm_df %>% filter(source == src & !is.na(score) & is.finite(score) & N_measured > 3)
    if (nrow(u_sub) == 0) next

    n_gse <- length(unique(u_sub$gse))
    if (n_gse < min_ds_req) next

    valid_comps <- u_sub$condition
    valid_gses <- u_sub$gse
    valid_ctrls <- u_sub$control
    valid_cov_cluster <- paste0(valid_gses, "_", ifelse(is.na(valid_ctrls) | valid_ctrls == "", "UnknownControl", valid_ctrls))

    # Exact ULM correlation -> Fisher-Z transformation
    df_u <- u_sub$N_measured - 2
    t_u <- u_sub$score
    r_u <- t_u / sqrt(t_u^2 + df_u)
    yi_u <- atanh(pmin(pmax(r_u, -0.999999), 0.999999))
    vi_u <- 1 / (u_sub$N_measured - 3)

    V_mat_u <- tryCatch(
      metafor::vcalc(vi = vi_u, cluster = valid_cov_cluster, obs = seq_along(yi_u), rho = 0.5),
      error = function(e) diag(vi_u, nrow = length(vi_u))
    )

    fit_u <- .fit_3lre_robust(yi = yi_u, V_mat = V_mat_u, gse_vec = valid_gses, comp_vec = valid_comps)
    fe_u <- tryCatch(metafor::rma.mv(yi = yi_u, V = V_mat_u), error = function(e) NULL)

    ulm_3lre_z <- NA_real_; ulm_3lre_se <- NA_real_; ulm_3lre_r <- NA_real_; ulm_3lre_p <- NA_real_
    ulm_i2_l3 <- NA_real_; ulm_i2_l2 <- NA_real_; ulm_tau2_l3 <- NA_real_; ulm_tau2_l2 <- NA_real_; ulm_het_p <- NA_real_
    if (!is.null(fit_u)) {
      s2 <- fit_u$re$sigma2
      ulm_tau2_l3 <- s2[1]; ulm_tau2_l2 <- s2[2]
      i2 <- .calc_pseudo_i2(vi_u, s2)
      ulm_i2_l3 <- i2[1]; ulm_i2_l2 <- i2[2]
      ulm_het_p <- fit_u$re$QEp
      est_obj <- if (!is.null(fit_u$rob)) fit_u$rob else fit_u$re
      ulm_3lre_z <- as.numeric(est_obj$beta)[1]
      ulm_3lre_se <- as.numeric(est_obj$se)[1]
      ulm_3lre_r <- tanh(ulm_3lre_z)
      ulm_3lre_p <- as.numeric(est_obj$pval)[1]
    }

    ulm_fe_z <- NA_real_; ulm_fe_se <- NA_real_; ulm_fe_r <- NA_real_; ulm_fe_p <- NA_real_
    if (!is.null(fe_u)) {
      ulm_fe_z <- as.numeric(fe_u$beta)[1]
      ulm_fe_se <- as.numeric(fe_u$se)[1]
      ulm_fe_r <- tanh(ulm_fe_z)
      ulm_fe_p <- as.numeric(fe_u$pval)[1]
    }

    med_targets <- median(vapply(valid_comps, function(cn) {
      val <- comp_target_counts[[cn]][src]
      if (is.na(val)) NA_real_ else as.numeric(val)
    }, numeric(1)), na.rm = TRUE)

    meta_rows[[i]] <- data.frame(
      source = src,
      condition = sample_label,
      n_comparisons = nrow(u_sub),
      n_datasets = n_gse,
      median_targets_mapped = med_targets,
      ULM_3LRE_FisherZ = ulm_3lre_z,
      ULM_3LRE_SE = ulm_3lre_se,
      ULM_3LRE_r = ulm_3lre_r,
      ULM_3LRE_pval = ulm_3lre_p,
      ULM_FE_FisherZ = ulm_fe_z,
      ULM_FE_SE = ulm_fe_se,
      ULM_FE_r = ulm_fe_r,
      ULM_FE_pval = ulm_fe_p,
      ULM_3LRE_I2_Level3 = ulm_i2_l3,
      ULM_3LRE_I2_Level2 = ulm_i2_l2,
      ULM_3LRE_tau2_Level3 = ulm_tau2_l3,
      ULM_3LRE_tau2_Level2 = ulm_tau2_l2,
      ULM_3LRE_Het_pval = ulm_het_p,
      n_pos_ULM = sum(u_sub$score > 0),
      n_neg_ULM = sum(u_sub$score < 0),
      n_sig_pos_ULM = sum(u_sub$score > 0 & u_sub$p_value < 0.05),
      n_sig_neg_ULM = sum(u_sub$score < 0 & u_sub$p_value < 0.05),
      stringsAsFactors = FALSE
    )
  }

  meta_df <- bind_rows(Filter(Negate(is.null), meta_rows))
  if (nrow(meta_df) == 0) {
    cat(sprintf("  ⚠ No iModulons met min_ds >= %d for meta-analysis pooling.\n", min_ds_req))
    return(NULL)
  }

  meta_df <- meta_df %>%
    mutate(
      ULM_3LRE_Z = ifelse(!is.na(ULM_3LRE_SE) & ULM_3LRE_SE > 0, ULM_3LRE_FisherZ / ULM_3LRE_SE, NA_real_),
      ULM_FE_Z = ifelse(!is.na(ULM_FE_SE) & ULM_FE_SE > 0, ULM_FE_FisherZ / ULM_FE_SE, NA_real_),
      ULM_3LRE_padj = p.adjust(ULM_3LRE_pval, method = "BH"),
      ULM_FE_padj = p.adjust(ULM_FE_pval, method = "BH")
    ) %>%
    left_join(imod_meta, by = "source") %>%
    arrange(ULM_3LRE_pval, desc(abs(ULM_3LRE_Z)))

  meta_csv <- file.path(out_dir, "iModulon_Meta_Analysis_Results.csv")
  write.csv(meta_df, meta_csv, row.names = FALSE)
  cat(sprintf("  ✓ Exported pooled iModulon meta-analysis (%d iModulons) -> %s\n", nrow(meta_df), meta_csv))

  # Barplot for pooled ULM meta-analysis (standardized Z = FisherZ / SE scale)
  .plot_imodulon_drivers(
    df = meta_df,
    score_col = "ULM_3LRE_Z",
    p_col = "ULM_3LRE_pval",
    padj_col = "ULM_3LRE_padj",
    title_str = paste("Top Meta-Analyzed iModulons (Per-Dataset ULM -> 3LRE):", sample_label),
    subtitle_str = "3-Level Random Effects on Fisher-Z (vcalc rho=0.5 + clubSandwich RVE)",
    x_lab = "Pooled ULM Activity (Z-score = Fisher Z / SE)",
    out_file = file.path(vis_dir, "Barplot_iModulon_Drivers_MetaAnalysis_ULM_3LRE.png")
  )

  # Heatmap of top iModulons across individual comparisons (ULM)
  top_src <- meta_df %>%
    arrange(ULM_3LRE_pval) %>%
    head(25) %>%
    pull(source)

  if (length(top_src) > 1) {
    heat_df <- all_per_df %>%
      filter(statistic == "ulm" & source %in% top_src) %>%
      mutate(
        source = factor(source, levels = rev(top_src)),
        score_capped = pmin(pmax(score, -4), 4)
      )

    p_heat <- ggplot(heat_df, aes(x = condition, y = source, fill = score_capped)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#1F78B4", mid = "white", high = "#D95F02", midpoint = 0, limits = c(-4, 4)) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 11)
      ) +
      labs(
        title = paste("Per-Dataset iModulon Activity Heatmap (Top 25, ULM):", sample_label),
        subtitle = "ULM t-scores capped at [-4, 4]; Biofilm vs Planktonic orientation",
        x = "Dataset Comparison",
        y = "iModulon / Regulator",
        fill = "ULM t"
      )
    heat_file <- file.path(vis_dir, "Heatmap_iModulon_PerDataset_Activities.png")
    ggsave(heat_file, plot = p_heat, width = max(10, length(per_comp_scores) * 0.28), height = 8, dpi = 300, limitsize = FALSE)
    cat("  ✓ Saved heatmap:", heat_file, "\n")
  }

  list(per_dataset = all_per_df, meta_summary = meta_df)
}

#' Compare Approach 1 (Direct Meta-Stats: ULM + Consensus) vs Approach 2 (Per-Dataset ULM -> Meta-Analysis)
#' @keywords internal
.compare_direct_vs_meta <- function(direct_results, meta_df, out_dir, run_consensus) {
  if (is.null(meta_df) || length(direct_results) == 0) return(invisible(NULL))

  cmp <- meta_df %>%
    select(
      source, imodulon_name, category, regulator,
      n_comparisons, n_datasets,
      MetaULM_3LRE_FisherZ = ULM_3LRE_FisherZ, MetaULM_3LRE_pval = ULM_3LRE_pval, MetaULM_3LRE_padj = ULM_3LRE_padj,
      MetaULM_FE_FisherZ = ULM_FE_FisherZ, MetaULM_FE_padj = ULM_FE_padj
    )

  stats_to_merge <- if (run_consensus) c("ulm", "consensus") else "ulm"
  for (m_lbl in names(direct_results)) {
    for (st in stats_to_merge) {
      d_sub <- direct_results[[m_lbl]] %>%
        filter(statistic == st) %>%
        select(source, score, p_value, padj)
      st_tag <- toupper(st)
      colnames(d_sub)[2:4] <- paste0("Direct_", m_lbl, "_", st_tag, c("_score", "_pval", "_padj"))
      cmp <- left_join(cmp, d_sub, by = "source")
    }
  }

  if ("Direct_MetaFE_ULM_padj" %in% colnames(cmp)) {
    sig_d <- !is.na(cmp$Direct_MetaFE_ULM_padj) & cmp$Direct_MetaFE_ULM_padj < 0.05
    sig_m <- !is.na(cmp$MetaULM_FE_padj) & cmp$MetaULM_FE_padj < 0.05
    same_dir <- sign(cmp$Direct_MetaFE_ULM_score) == sign(cmp$MetaULM_FE_FisherZ)
    cmp$Concordance_DirectFE_ULM_vs_MetaULM_FE <- ifelse(
      sig_d & sig_m & same_dir,
      ifelse(cmp$Direct_MetaFE_ULM_score > 0, "Both_Activated", "Both_Repressed"),
      ifelse(sig_d & sig_m & !same_dir, "Discordant",
        ifelse(sig_d, "Direct_FE_ULM_only", ifelse(sig_m, "Meta_ULM_FE_only", "Neither")))
    )
  }

  cmp_file <- file.path(out_dir, "iModulon_Direct_vs_Meta_Comparison.csv")
  write.csv(cmp, cmp_file, row.names = FALSE)
  cat("  ✓ Exported cross-approach comparison ->", cmp_file, "\n")
  cmp
}

#' Run iModulon Activity Inference using decoupleR (Direct + Per-Dataset Meta-Analysis)
#'
#' @param res_df Data frame of meta-analysis results (res_df_export from DESeq2 script).
#' @param lfc_list Optional named list of per-comparison LFC vectors (for Approach 2).
#' @param se_list Optional named list of per-comparison SE vectors (for Approach 2).
#' @param comparison_gse_map Optional named character vector mapping comparison -> GSE/study cluster.
#' @param comparison_control_map Optional named character vector mapping comparison -> control group.
#' @param neff_list Optional named list of per-comparison effective sample sizes.
#' @param min_ds Minimum GSE datasets required to meta-analyze an iModulon (default: 2).
#' @param sample_label Condition label or sample_dir name (string after "sample_info-" is used).
#' @param m_file Path to iModulonDB M matrix CSV.
#' @param table_file Path to iModulonDB iModulon metadata CSV.
#' @param thresholds_file Path to per-iModulon thresholds CSV.
#' @param gene_info_file Path to staph_precise165-gene_info.csv.
#' @param trn_file Path to staph_precise165-trn.csv.
#' @param orthologue_file Path to OrthologueTable.xlsx.
#' @param annote_dir Directory containing GeneSpecificInformation_*.xlsx files.
#' @param operon_file Path to operon_gene_mapping.csv.
#' @param gaf_file Path to QuickGO GAF file.
#' @param out_dir Output directory for CSV tables.
#' @param vis_dir Output directory for visualization plots.
#' @param weight_threshold Fallback absolute cutoff for M matrix weights (default: 0.08).
#' @param min_targets Minimum target genes required to evaluate an iModulon (default: 5).
#' @param run_consensus Logical; compute consensus across ULM, WMEAN, and MLM (default: TRUE).
run_imodulon_analysis <- function(res_df,
                                  lfc_list = NULL,
                                  se_list = NULL,
                                  comparison_gse_map = NULL,
                                  comparison_control_map = NULL,
                                  neff_list = NULL,
                                  min_ds = 2,
                                  sample_label = "Biofilm_vs_Planktonic",
                                  m_file = "imodulon_data/staph_precise165-M.csv",
                                  table_file = "imodulon_data/staph_precise165-im_table.csv",
                                  thresholds_file = "imodulon_data/staph_precise165-M_thresholds.csv",
                                  gene_info_file = "imodulon_data/staph_precise165-gene_info.csv",
                                  trn_file = "imodulon_data/staph_precise165-trn.csv",
                                  orthologue_file = "OrthologueTable.xlsx",
                                  annote_dir = "annote",
                                  operon_file = "operon_gene_mapping.csv",
                                  gaf_file = "QuickGO-annotations-1764241851574-20251127.gaf",
                                  out_dir = "results",
                                  vis_dir = file.path(out_dir, "Visualizations_Independent"),
                                  weight_threshold = 0.08,
                                  min_targets = 5,
                                  run_consensus = TRUE) {
  cat("\n=======================================================\n")
  cat("RUNNING IMODULON / DECOUPLER REGULATORY META-ANALYSIS\n")
  cat("=======================================================\n")

  if (!file.exists(m_file)) stop(paste("CRITICAL ERROR: iModulon M matrix file not found at:", m_file))
  if (!file.exists(table_file)) stop(paste("CRITICAL ERROR: iModulon table metadata file not found at:", table_file))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!dir.exists(vis_dir)) dir.create(vis_dir, recursive = TRUE)

  # Normalize condition label to the string after "sample_info-"
  clean_label <- sub("^sample_info-", "", basename(as.character(sample_label)))
  if (is.na(clean_label) || clean_label == "") clean_label <- "Biofilm_vs_Planktonic"
  cat("  Condition column label:", clean_label, "\n")

  # 1. Curate iModulon network
  curated <- .curate_imodulon_network(m_file, table_file, thresholds_file, weight_threshold)

  # 2. Harmonize cluster_ids <-> M-matrix targets & report mapping %
  if (!"cluster_id" %in% colnames(res_df)) res_df$cluster_id <- rownames(res_df)
  cluster_df <- res_df %>% select(any_of(c("cluster_id", "Representative_Gene", "All_Genes_Locus_Tags")))
  mapped <- .build_mapped_network(
    net_raw = curated$net_raw,
    cluster_df = cluster_df,
    m_gene_ids = curated$m_gene_ids,
    gene_info_file = gene_info_file,
    trn_file = trn_file,
    orthologue_file = orthologue_file,
    annote_dir = annote_dir,
    operon_file = operon_file,
    gaf_file = gaf_file,
    out_dir = out_dir
  )

  # 3. Approach 1: Direct decoupleR on pooled meta-analysis statistics
  direct_res <- .run_approach1_direct(
    res_df = res_df,
    net_mapped = mapped$net_mapped,
    imod_meta = curated$imod_meta,
    sample_label = clean_label,
    out_dir = out_dir,
    vis_dir = vis_dir,
    min_targets = min_targets,
    run_consensus = run_consensus
  )

  # 4. Approach 2: Per-dataset ULM -> 3LRE / FE meta-analysis
  meta_res <- .run_approach2_per_dataset_meta(
    lfc_list = lfc_list,
    se_list = se_list,
    comparison_gse_map = comparison_gse_map,
    comparison_control_map = comparison_control_map,
    neff_list = neff_list,
    min_ds = min_ds,
    net_mapped = mapped$net_mapped,
    imod_meta = curated$imod_meta,
    sample_label = clean_label,
    out_dir = out_dir,
    vis_dir = vis_dir,
    min_targets = min_targets
  )

  # 5. Cross-approach comparison
  cmp_res <- .compare_direct_vs_meta(direct_res, meta_res$meta_summary, out_dir, run_consensus)

  cat("\n✓ iModulon decoupleR analysis finished successfully.\n")
  list(
    direct = direct_res,
    per_dataset = meta_res$per_dataset,
    meta = meta_res$meta_summary,
    comparison = cmp_res,
    mapping_stats = mapped$mapping_stats
  )
}

#' Standalone runner on an existing results_sample_info-* directory
run_imodulon_from_results_dir <- function(results_dir = "results_sample_info-all", min_ds = 2, min_targets = 5) {
  res_file <- file.path(results_dir, "DESeq2_contrast_results.csv")
  coef_file <- file.path(results_dir, "dataset_coefficients.csv")
  if (!file.exists(res_file)) stop("Missing DESeq2_contrast_results.csv in: ", results_dir)

  res_df <- read.csv(res_file, stringsAsFactors = FALSE)
  sample_label <- sub("^results_(sample_info-)?", "", basename(results_dir))

  lfc_list <- NULL; se_list <- NULL; comp_gse <- NULL; comp_ctrl <- NULL
  if (file.exists(coef_file)) {
    coef_df <- read.csv(coef_file, check.names = FALSE, stringsAsFactors = FALSE)
    lfc_cols <- grep("_LFC_adj$", colnames(coef_df), value = TRUE)
    comp_names <- sub("_LFC_adj$", "", lfc_cols)

    lfc_list <- list(); se_list <- list()
    comp_gse <- character(0); comp_ctrl <- character(0)
    for (cn in comp_names) {
      se_col <- paste0(cn, "_SE")
      if (!se_col %in% colnames(coef_df)) next
      lfc_list[[cn]] <- setNames(coef_df[[paste0(cn, "_LFC_adj")]], coef_df$cluster_id)
      se_list[[cn]] <- setNames(coef_df[[se_col]], coef_df$cluster_id)
      comp_gse[cn] <- sub("_.*", "", cn)
      comp_ctrl[cn] <- if (grepl("_vs_", cn)) sub(".*_vs_", "", cn) else NA_character_
    }
  }

  run_imodulon_analysis(
    res_df = res_df,
    lfc_list = lfc_list,
    se_list = se_list,
    comparison_gse_map = comp_gse,
    comparison_control_map = comp_ctrl,
    min_ds = min_ds,
    sample_label = sample_label,
    out_dir = results_dir,
    vis_dir = file.path(results_dir, "Visualizations_Independent"),
    min_targets = min_targets,
    run_consensus = TRUE
  )
}

# Allow direct CLI execution for standalone testing on existing results directories
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  target_dir <- if (length(args) >= 1) args[1] else "results_sample_info-all"
  min_ds_arg <- if (length(args) >= 2) as.integer(args[2]) else 12L
  min_t_arg <- if (length(args) >= 3) as.integer(args[3]) else 5L
  run_imodulon_from_results_dir(target_dir, min_ds = min_ds_arg, min_targets = min_t_arg)
}
