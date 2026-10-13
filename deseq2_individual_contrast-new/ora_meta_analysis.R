run_ora_meta_analysis <- function(gaf_file, lfc_list, pval_list, se_list, cts, cluster_annotation, comparison_gse_map, comparison_control_map = NULL, min_ds, out_dir, ora_padj = 0.05, ora_lfc = 1.0) {
  cat("\n=== RUNNING OVER-REPRESENTATION META-ANALYSIS ===\n")

  # Convert a two-sided p-value to a one-sided p-value testing estimate > 0.
  # Positive estimate: halve p. Negative estimate: 1 - p/2 (never significant).
  one_sided_p <- function(est, p) ifelse(est > 0, p / 2, 1 - p / 2)
  if (!is.null(gaf_file) && file.exists(gaf_file)) {
    cat("Loading GAF file:", gaf_file, "\n")

    # --- 1. Setup and Package Loading ---
    # Ensure required Bioconductor and CRAN packages are installed and loaded
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
    }
    if (!requireNamespace("GO.db", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("GO.db", update = FALSE, ask = FALSE)
    }
    if (!requireNamespace("clubSandwich", quietly = TRUE)) {
      install.packages("clubSandwich", repos = "https://cloud.r-project.org", quiet = TRUE)
    }

    suppressPackageStartupMessages(library(clusterProfiler))
    suppressPackageStartupMessages(library(GO.db))
    suppressPackageStartupMessages(library(clubSandwich))

    # --- 2. GAF Parsing and Cluster ID Mapping ---
    # Read the custom GAF file and map dataset cluster IDs to standard gene symbols
    gaf <- read.delim(gaf_file, comment.char = "!", header = FALSE, stringsAsFactors = FALSE, quote = "")
    term2gene <- unique(data.frame(term = gaf$V5, gene = gaf$V3, stringsAsFactors = FALSE))

    # The GAF contains only direct annotations. GO is a DAG, so a gene annotated to a
    # term is implicitly annotated to all its ancestors. Propagate them with buildGOmap().
    n_direct <- nrow(term2gene)
    term2gene <- clusterProfiler::buildGOmap(term2gene)
    term2gene <- unique(data.frame(term = as.character(term2gene[[1]]), gene = as.character(term2gene[[2]]), stringsAsFactors = FALSE))
    # Drop the ontology root terms (BP/MF/CC), which contain every gene and carry no information
    term2gene <- term2gene[grepl("^GO:", term2gene$term) & !term2gene$term %in% c("GO:0008150", "GO:0003674", "GO:0005575"), ]
    cat(sprintf("  Propagated GO annotations: %d direct -> %d direct+indirect term-gene pairs\n", n_direct, nrow(term2gene)))

    annotated_genes <- unique(gaf$V3)

    cluster_id_to_symbol <- sapply(cluster_annotation$All_Genes_Locus_Tags, function(tags) {
      if (is.na(tags) || tags == "") {
        return(NA_character_)
      }
      tag_list <- strsplit(as.character(tags), "/")[[1]]
      matches <- tag_list[tag_list %in% annotated_genes]
      if (length(matches) > 0) {
        return(matches[1])
      } else {
        return(NA_character_)
      }
    })
    names(cluster_id_to_symbol) <- cluster_annotation$cluster_id

    goterms <- as.data.frame(GO.db::GOTERM)
    term2name <- unique(data.frame(term = goterms$go_id, name = goterms$Term, stringsAsFactors = FALSE))
    term2ontology <- unique(data.frame(term = goterms$go_id, ontology = goterms$Ontology, stringsAsFactors = FALSE))

    ora_logor_list <- list()
    ora_se_list <- list()
    ora_z_list <- list()
    ora_z_se_list <- list()
    gsea_z_list <- list()
    ora_comparisons <- c()
    gsea_comparisons <- c()
    ora_gse_map <- c()

    # --- 3. Compute ORA for Each Comparison ---
    # Iterate through each dataset contrast to perform dataset-specific over-representation analysis
    for (comp_name in names(lfc_list)) {
      gse_id <- comparison_gse_map[comp_name]

      pvals <- pval_list[[comp_name]]
      lfcs <- lfc_list[[comp_name]]
      ses <- se_list[[comp_name]]

      valid_mask <- !is.na(pvals) & !is.na(lfcs)
      universe_clusters <- rownames(cts)[valid_mask]

      if (length(universe_clusters) == 0) next

      dataset_padj <- p.adjust(pvals[valid_mask], method = "BH")
      de_mask <- dataset_padj < ora_padj & abs(lfcs[valid_mask]) >= ora_lfc
      de_clusters <- universe_clusters[de_mask]

      universe_symbols <- unique(cluster_id_to_symbol[universe_clusters])
      universe_symbols <- universe_symbols[!is.na(universe_symbols) & universe_symbols != ""]

      de_symbols <- unique(cluster_id_to_symbol[de_clusters])
      de_symbols <- de_symbols[!is.na(de_symbols) & de_symbols != ""]

      cat(sprintf(
        "  [%s] Eligible genes: %d | Mapped universe genes: %d | Mapped DE genes for ORA: %d\n",
        comp_name, length(universe_clusters), length(universe_symbols), length(de_symbols)
      ))

      if (length(de_symbols) == 0) {
        cat("  - Skipping ORA for", comp_name, "due to 0 mapped DE genes.\n")
        next
      }

      ora_res <- tryCatch(
        {
          enricher(
            gene = de_symbols,
            universe = universe_symbols,
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            minGSSize = 5,
            maxGSSize = 500,
            TERM2GENE = term2gene,
            TERM2NAME = term2name
          )
        },
        error = function(e) NULL
      )

      # --- 4. Calculate Log Odds Ratios ---
      # If ORA returns results, compute the continuous Log Odds Ratio (logOR) and Standard Error (SE)
      # using the Haldane-Anscombe correction (+0.5) to avoid zeroes
      if (!is.null(ora_res) && nrow(ora_res@result) > 0) {
        res_df <- ora_res@result

        parse_num <- function(x) as.numeric(sapply(strsplit(as.character(x), "/"), `[`, 1))
        parse_den <- function(x) as.numeric(sapply(strsplit(as.character(x), "/"), `[`, 2))

        a <- parse_num(res_df$GeneRatio)
        de_size <- parse_den(res_df$GeneRatio)
        b <- de_size - a

        bg_num <- parse_num(res_df$BgRatio)
        bg_den <- parse_den(res_df$BgRatio)
        c_val <- bg_num - a
        d <- bg_den - de_size - c_val

        a_adj <- a + 0.5
        b_adj <- b + 0.5
        c_adj <- c_val + 0.5
        d_adj <- d + 0.5

        logOR <- log((a_adj * d_adj) / (b_adj * c_adj))
        se_logOR <- sqrt(1 / a_adj + 1 / b_adj + 1 / c_adj + 1 / d_adj)

        z_scores <- qnorm(res_df$pvalue, lower.tail = FALSE)
        z_scores[is.infinite(z_scores)] <- 8.2
        se_z <- rep(1, length(z_scores))

        ora_logor_list[[comp_name]] <- setNames(logOR, res_df$ID)
        ora_se_list[[comp_name]] <- setNames(se_logOR, res_df$ID)
        ora_z_list[[comp_name]] <- setNames(z_scores, res_df$ID)
        ora_z_se_list[[comp_name]] <- setNames(se_z, res_df$ID)
        ora_comparisons <- c(ora_comparisons, comp_name)
        ora_gse_map[comp_name] <- gse_id
      }

      # --- 4. Compute GSEA ---
      valid_rank_mask <- !is.na(lfcs) & !is.na(ses) & ses > 0
      rank_clusters <- universe_clusters[valid_rank_mask]
      rank_symbols <- cluster_id_to_symbol[rank_clusters]

      t_stats <- lfcs[valid_rank_mask] / ses[valid_rank_mask]
      names(t_stats) <- rank_symbols

      valid_sym_mask <- !is.na(names(t_stats)) & names(t_stats) != "" & is.finite(t_stats)
      t_stats <- t_stats[valid_sym_mask]

      if (length(t_stats) > 0) {
        if (any(duplicated(names(t_stats)))) {
          t_stats <- tapply(t_stats, names(t_stats), mean, na.rm = TRUE)
          t_stats <- setNames(as.numeric(t_stats), names(t_stats))
        }
        t_stats <- sort(t_stats, decreasing = TRUE)

        gsea_res <- tryCatch(
          {
            clusterProfiler::GSEA(
              geneList = t_stats,
              pvalueCutoff = 1,
              minGSSize = 5,
              maxGSSize = 500,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              verbose = FALSE,
              eps = 0
            )
          },
          error = function(e) {
            cat("  - GSEA failed for", comp_name, ":", e$message, "\n")
            NULL
          }
        )

        if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
          res_df_g <- gsea_res@result
          safe_p <- pmax(res_df_g$pvalue, 1e-10)
          z_dir <- qnorm(safe_p, lower.tail = FALSE)
          gsea_z_list[[comp_name]] <- setNames(z_dir, res_df_g$ID)
          gsea_comparisons <- c(gsea_comparisons, comp_name)
        }
      }
    }

    # --- 5. Aggregate Effect Sizes ---
    # Combine the logOR and SE from all individual comparisons into master matrices
    if (length(ora_logor_list) > 0) {
      all_go_ids <- unique(unlist(lapply(ora_logor_list, names)))

      ora_yi <- matrix(NA, nrow = length(all_go_ids), ncol = length(ora_logor_list), dimnames = list(all_go_ids, ora_comparisons))
      ora_vi <- matrix(NA, nrow = length(all_go_ids), ncol = length(ora_se_list), dimnames = list(all_go_ids, ora_comparisons))
      ora_z_yi <- matrix(NA, nrow = length(all_go_ids), ncol = length(ora_z_list), dimnames = list(all_go_ids, ora_comparisons))
      ora_z_vi <- matrix(NA, nrow = length(all_go_ids), ncol = length(ora_z_se_list), dimnames = list(all_go_ids, ora_comparisons))

      for (comp_name in ora_comparisons) {
        ora_yi[names(ora_logor_list[[comp_name]]), comp_name] <- ora_logor_list[[comp_name]]
        ora_vi[names(ora_se_list[[comp_name]]), comp_name] <- ora_se_list[[comp_name]]^2
        ora_z_yi[names(ora_z_list[[comp_name]]), comp_name] <- ora_z_list[[comp_name]]
        ora_z_vi[names(ora_z_se_list[[comp_name]]), comp_name] <- ora_z_se_list[[comp_name]]^2
      }

      cat("  Meta-analyzing", length(all_go_ids), "GO terms across", length(ora_comparisons), "comparisons...\n")

      ora_meta_res <- list()
      debug_printed <- FALSE
      for (i in seq_along(all_go_ids)) {
        go_id <- all_go_ids[i]
        yi_go <- ora_yi[i, ]
        vi_go <- ora_vi[i, ]
        z_yi_go <- ora_z_yi[i, ]
        z_vi_go <- ora_z_vi[i, ]

        # --- 6. Meta-Analyze Each GO Term ---
        # Fit 3-Level Random Effects (w/ clubSandwich RVE) and Fixed Effects models
        valid_go <- !is.na(yi_go) & !is.na(vi_go)
        valid_gse_go <- ora_gse_map[valid_go]
        valid_comp_go <- ora_comparisons[valid_go]

        if (length(unique(valid_gse_go)) >= max(2, min_ds)) {
          # Determine clustering for shared controls to assign rho = 0.5
          if (!is.null(comparison_control_map)) {
            valid_control_go <- comparison_control_map[valid_comp_go]
            # Use GSE + Control as the cluster. If control is NA, group them together in an 'UnknownControl' cluster for 0.5 correlation
            valid_cov_cluster <- paste0(valid_gse_go, "_", ifelse(is.na(valid_control_go), "UnknownControl", valid_control_go))
          } else {
            # If no map is provided, assume rho=0.5 between comparisons in the same GSE dataset
            valid_cov_cluster <- valid_gse_go
          }

          # Impute covariance (rho = 0.5) for comparisons sharing the same dataset/control
          V_mat_logor <- tryCatch(
            metafor::vcalc(vi = vi_go[valid_go], cluster = valid_cov_cluster, obs = 1:sum(valid_go), rho = 0.5),
            error = function(e) vi_go[valid_go]
          )

          # Z-score Meta-Analysis Covariance Matrix
          V_mat_z <- tryCatch(
            metafor::vcalc(vi = z_vi_go[valid_go], cluster = valid_cov_cluster, obs = 1:sum(valid_go), rho = 0.5),
            error = function(e) z_vi_go[valid_go]
          )

          if (!debug_printed) {
            cat("\n--- DEBUG: Covariance Matrices for GO term", go_id, "---\n")
            cat("Comparisons:\n")
            print(valid_comp_go)
            cat("\nClusters (valid_cov_cluster):\n")
            print(valid_cov_cluster)
            cat("\nCorrelation Matrix (derived from V_mat_logor):\n")
            print(cov2cor(V_mat_logor))
            cat("\nCorrelation Matrix (derived from V_mat_z):\n")
            print(cov2cor(V_mat_z))
            cat("-----------------------------------------------\n\n")
            debug_printed <- TRUE
          }

          # LogOR Meta-Analysis
          res_re <- tryCatch(
            metafor::rma.mv(yi = yi_go[valid_go], V = V_mat_logor, random = ~ 1 | valid_gse_go / valid_comp_go, method = "REML"),
            error = function(e) {
              cat(sprintf("  - [logOR] 3-level RE model failed for %s: %s\n", go_id, e$message))
              NULL
            }
          )

          if (!is.null(res_re)) {
            res_re <- tryCatch(
              metafor::robust(res_re, cluster = valid_gse_go, clubSandwich = TRUE),
              error = function(e) {
                cat(sprintf("  - [logOR] clubSandwich RVE failed for %s: %s\n", go_id, e$message))
                NULL
              }
            )
          }

          res_fe <- tryCatch(
            metafor::rma.mv(yi = yi_go[valid_go], V = V_mat_logor),
            error = function(e) {
              cat(sprintf("  - [logOR FE] model failed for %s: %s\n", go_id, e$message))
              NULL
            }
          )

          # Z-score Meta-Analysis Models

          res_z_re <- tryCatch(
            metafor::rma.mv(yi = z_yi_go[valid_go], V = V_mat_z, random = ~ 1 | valid_gse_go / valid_comp_go, method = "REML"),
            error = function(e) {
              cat(sprintf("  - [Z-score] 3-level RE model failed for %s: %s\n", go_id, e$message))
              NULL
            }
          )

          if (!is.null(res_z_re)) {
            res_z_re <- tryCatch(
              metafor::robust(res_z_re, cluster = valid_gse_go, clubSandwich = TRUE),
              error = function(e) {
                cat(sprintf("  - [Z-score] clubSandwich RVE failed for %s: %s\n", go_id, e$message))
                NULL
              }
            )
          }

          res_z_fe <- tryCatch(
            metafor::rma.mv(yi = z_yi_go[valid_go], V = V_mat_z),
            error = function(e) {
              cat(sprintf("  - [Z-score FE] model failed for %s: %s\n", go_id, e$message))
              NULL
            }
          )

          if (!is.null(res_re) && !is.null(res_fe) && !is.null(res_z_re) && !is.null(res_z_fe)) {
            ora_meta_res[[go_id]] <- data.frame(
              ID = go_id,
              Ontology = term2ontology$ontology[match(go_id, term2ontology$term)],
              Description = term2name$name[match(go_id, term2name$term)],
              logOR_RE_estimate = res_re$beta[1],
              logOR_RE_SE = res_re$se[1],
              logOR_RE_pval = one_sided_p(res_re$beta[1], res_re$pval[1]),
              logOR_FE_estimate = res_fe$beta[1],
              logOR_FE_SE = res_fe$se[1],
              logOR_FE_pval = one_sided_p(res_fe$beta[1], res_fe$pval[1]),
              Z_RE_estimate = res_z_re$beta[1],
              Z_RE_SE = res_z_re$se[1],
              Z_RE_pval = one_sided_p(res_z_re$beta[1], res_z_re$pval[1]),
              Z_FE_estimate = res_z_fe$beta[1],
              Z_FE_SE = res_z_fe$se[1],
              Z_FE_pval = one_sided_p(res_z_fe$beta[1], res_z_fe$pval[1]),
              N_comparisons = sum(valid_go),
              N_datasets = length(unique(valid_gse_go)),
              stringsAsFactors = FALSE
            )
          }
        }
      }

      # --- 7. Save and Export Results ---
      # Combine all meta-analyzed GO terms, sort by Random-Effects p-value, and export to CSV
      if (length(ora_meta_res) > 0) {
        ora_meta_df <- do.call(rbind, ora_meta_res)

        # Calculate FDR adjusted p-values for the meta-analysis results
        ora_meta_df$logOR_RE_padj <- p.adjust(ora_meta_df$logOR_RE_pval, method = "BH")
        ora_meta_df$logOR_FE_padj <- p.adjust(ora_meta_df$logOR_FE_pval, method = "BH")
        ora_meta_df$Z_RE_padj <- p.adjust(ora_meta_df$Z_RE_pval, method = "BH")
        ora_meta_df$Z_FE_padj <- p.adjust(ora_meta_df$Z_FE_pval, method = "BH")

        # Reorder columns to put padj next to pval
        ora_meta_df <- ora_meta_df[, c(
          "ID", "Ontology", "Description",
          "logOR_RE_estimate", "logOR_RE_SE", "logOR_RE_pval", "logOR_RE_padj",
          "logOR_FE_estimate", "logOR_FE_SE", "logOR_FE_pval", "logOR_FE_padj",
          "Z_RE_estimate", "Z_RE_SE", "Z_RE_pval", "Z_RE_padj",
          "Z_FE_estimate", "Z_FE_SE", "Z_FE_pval", "Z_FE_padj",
          "N_comparisons", "N_datasets"
        )]

        ora_meta_df <- ora_meta_df[order(ora_meta_df$logOR_RE_pval), ]
        write.csv(ora_meta_df, file.path(out_dir, "ORA_Meta_Analysis_Results.csv"), row.names = FALSE)
        cat("  - Saved ORA Meta-Analysis Results to", file.path(out_dir, "ORA_Meta_Analysis_Results.csv"), "\n")
      } else {
        cat("  - Not enough valid datasets/comparisons to run meta-analysis on any GO term.\n")
      }
    } else {
      cat("  - No ORA results could be computed for any comparison.\n")
    }

    # --- 8. Aggregate and Meta-Analyze GSEA ---
    if (length(gsea_z_list) > 0) {
      all_gsea_ids <- unique(unlist(lapply(gsea_z_list, names)))
      gsea_yi <- matrix(NA, nrow = length(all_gsea_ids), ncol = length(gsea_z_list), dimnames = list(all_gsea_ids, gsea_comparisons))
      for (comp_name in gsea_comparisons) {
        gsea_yi[names(gsea_z_list[[comp_name]]), comp_name] <- gsea_z_list[[comp_name]]
      }

      cat("  Meta-analyzing", length(all_gsea_ids), "GO terms (GSEA) across", length(gsea_comparisons), "comparisons...\n")
      gsea_meta_res <- list()
      for (i in seq_along(all_gsea_ids)) {
        go_id <- all_gsea_ids[i]
        yi_go <- gsea_yi[i, ]
        valid_go <- !is.na(yi_go)
        valid_gse_go <- ora_gse_map[valid_go]
        valid_comp_go <- gsea_comparisons[valid_go]

        if (length(unique(valid_gse_go)) >= max(2, min_ds)) {
          if (!is.null(comparison_control_map)) {
            valid_control_go <- comparison_control_map[valid_comp_go]
            valid_cov_cluster <- paste0(valid_gse_go, "_", ifelse(is.na(valid_control_go), "UnknownControl", valid_control_go))
          } else {
            # If no map is provided, assume rho=0.5 between comparisons in the same GSE dataset
            valid_cov_cluster <- valid_gse_go
          }

          V_mat_z <- tryCatch(
            metafor::vcalc(vi = rep(1, sum(valid_go)), cluster = valid_cov_cluster, obs = 1:sum(valid_go), rho = 0.5),
            error = function(e) rep(1, sum(valid_go))
          )

          res_re <- tryCatch(
            metafor::rma.mv(yi = yi_go[valid_go], V = V_mat_z, random = ~ 1 | valid_gse_go / valid_comp_go, method = "REML"),
            error = function(e) NULL
          )
          if (!is.null(res_re)) {
            res_re <- tryCatch(metafor::robust(res_re, cluster = valid_gse_go, clubSandwich = TRUE), error = function(e) NULL)
          }

          res_fe <- tryCatch(
            metafor::rma.mv(yi = yi_go[valid_go], V = V_mat_z),
            error = function(e) NULL
          )

          if (!is.null(res_re) && !is.null(res_fe)) {
            gsea_meta_res[[go_id]] <- data.frame(
              ID = go_id,
              Ontology = term2ontology$ontology[match(go_id, term2ontology$term)],
              Description = term2name$name[match(go_id, term2name$term)],
              GSEA_Z_RE_estimate = res_re$beta[1],
              GSEA_Z_RE_SE = res_re$se[1],
              GSEA_Z_RE_pval = one_sided_p(res_re$beta[1], res_re$pval[1]),
              GSEA_Z_FE_estimate = res_fe$beta[1],
              GSEA_Z_FE_SE = res_fe$se[1],
              GSEA_Z_FE_pval = one_sided_p(res_fe$beta[1], res_fe$pval[1]),
              N_comparisons = sum(valid_go),
              N_datasets = length(unique(valid_gse_go)),
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(gsea_meta_res) > 0) {
        gsea_df <- do.call(rbind, gsea_meta_res)
        gsea_df$GSEA_Z_RE_padj <- p.adjust(gsea_df$GSEA_Z_RE_pval, method = "BH")
        gsea_df$GSEA_Z_FE_padj <- p.adjust(gsea_df$GSEA_Z_FE_pval, method = "BH")
        gsea_df <- gsea_df[order(gsea_df$GSEA_Z_RE_pval), ]
        write.csv(gsea_df, file.path(out_dir, "GSEA_Meta_Analysis_Results.csv"), row.names = FALSE)
        cat("  - Saved GSEA Meta-Analysis Results to", file.path(out_dir, "GSEA_Meta_Analysis_Results.csv"), "\n")
      }
    } else {
      cat("  - No GSEA results could be computed for any comparison.\n")
    }
  } else {
    cat("  - GAF file not found or not provided. Skipping ORA meta-analysis.\n")
  }
}
