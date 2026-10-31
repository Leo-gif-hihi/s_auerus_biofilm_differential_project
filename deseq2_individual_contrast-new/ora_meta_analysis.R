# Reporting-only redundancy reduction. Called AFTER BH adjustment; it never alters the
# tested terms or their p-values/padj. Significant terms (padj < sig_padj) are clustered
# per ontology by GOSemSim (Wang graph-based similarity, hierarchical average linkage,
# cut at similarity >= sim_cutoff). The most significant term (lowest padj, ties broken
# by largest estimate) represents each cluster.
cluster_significant_go_terms <- function(res_df, sig_padj = 0.05, sim_cutoff = 0.7, padj_col = "logOR_RE_padj", est_col = "logOR_RE_estimate", measure = "Wang") {
  if (!requireNamespace("GOSemSim", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("GOSemSim", update = FALSE, ask = FALSE)
  }
  sig <- res_df[!is.na(res_df[[padj_col]]) & res_df[[padj_col]] < sig_padj, , drop = FALSE]
  if (nrow(sig) == 0) return(NULL)

  out <- list()
  for (ont in unique(sig$Ontology)) {
    sub <- sig[sig$Ontology == ont, , drop = FALSE]
    sub <- sub[order(sub[[padj_col]], -sub[[est_col]]), , drop = FALSE]
    if (nrow(sub) == 1) {
      sub$Cluster <- paste0(ont, "_1")
      sub$Cluster_Size <- 1L
      sub$Is_Representative <- TRUE
      out[[ont]] <- sub
      next
    }
    sem <- GOSemSim::godata(OrgDb = NULL, ont = ont, computeIC = FALSE)
    sim <- suppressMessages(GOSemSim::mgoSim(sub$ID, sub$ID, semData = sem, measure = measure, combine = NULL))
    sim[is.na(sim)] <- 0 # obsolete/unmappable terms -> own singleton cluster
    diag(sim) <- 1
    hc <- hclust(as.dist(1 - sim), method = "average")
    cl <- cutree(hc, h = 1 - sim_cutoff)
    sub$Cluster <- paste0(ont, "_", unname(cl[match(sub$ID, rownames(sim))]))
    sub$Cluster_Size <- as.integer(table(sub$Cluster)[sub$Cluster])
    # sub is sorted by significance, so the first term in each cluster is the representative
    sub$Is_Representative <- !duplicated(sub$Cluster)
    out[[ont]] <- sub
  }
  do.call(rbind, out)
}

run_ora_meta_analysis <- function(gaf_file, lfc_list, pval_list, se_list, cts, cluster_annotation, comparison_gse_map, comparison_control_map = NULL, min_ds, out_dir, ora_padj = 0.05, ora_lfc = 1.0, exclude_evidence = c("ND"), min_gs_size = 5, max_gs_size = 500, directions = c("all", "up", "down"), sem_sim_cutoff = 0.7, report_padj = 0.05) {
  cat("\n=== RUNNING OVER-REPRESENTATION META-ANALYSIS ===\n")

  # Validate directions
  directions <- match.arg(directions, c("all", "up", "down"), several.ok = TRUE)
  directions <- unique(directions)
  cat(sprintf("Directions to analyze: %s\n", paste(directions, collapse = ", ")))

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

    # Filter GAF annotations (Issue 5):
    # - Exclude negative qualifiers in V4 (e.g. 'NOT', 'NOT|enables', 'NOT|involved_in')
    # - Exclude uninformative evidence codes in V7 (e.g. 'ND' = No biological Data available)
    n_raw <- nrow(gaf)
    gaf <- gaf[!grepl("(^|\\|)NOT(\\||$)", gaf$V4), ]
    if (!is.null(exclude_evidence) && length(exclude_evidence) > 0) {
      gaf <- gaf[!gaf$V7 %in% exclude_evidence, ]
    }
    cat(sprintf("  GAF filtering: %d raw -> %d annotations retained (filtered NOT qualifiers; excluded evidence: %s)\n",
                n_raw, nrow(gaf), if (length(exclude_evidence) > 0) paste(exclude_evidence, collapse = ", ") else "none"))

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

    # --- Global Gene-Set Size Filtering (Issue 2) ---
    # Apply min_gs_size and max_gs_size globally against the study universe
    # so that the set of tested GO terms is uniform across all comparisons and does not
    # vary dynamically depending on individual comparison universes.
    all_genes <- rownames(cts)
    global_symbols <- unique(cluster_id_to_symbol[all_genes])
    global_symbols <- global_symbols[!is.na(global_symbols) & global_symbols != ""]

    term2gene_study <- term2gene[term2gene$gene %in% global_symbols, ]
    term_sizes <- table(term2gene_study$term)
    valid_terms <- names(term_sizes)[term_sizes >= min_gs_size & term_sizes <= max_gs_size]
    term2gene <- term2gene_study[term2gene_study$term %in% valid_terms, ]
    cat(sprintf("  Global gene-set size filtering (%d-%d): %d terms retained across %d study genes\n",
                min_gs_size, max_gs_size, length(valid_terms), length(global_symbols)))

    goterms <- as.data.frame(GO.db::GOTERM)
    term2name <- unique(data.frame(term = goterms$go_id, name = goterms$Term, stringsAsFactors = FALSE))
    term2ontology <- unique(data.frame(term = goterms$go_id, ontology = goterms$Ontology, stringsAsFactors = FALSE))

    ora_logor_list <- setNames(lapply(directions, function(x) list()), directions)
    ora_se_list <- setNames(lapply(directions, function(x) list()), directions)
    ora_comparisons <- setNames(lapply(directions, function(x) character(0)), directions)
    ora_gse_map <- setNames(lapply(directions, function(x) character(0)), directions)
    ora_contributing_genes <- setNames(lapply(directions, function(x) list()), directions)

    abs_ora_lfc <- abs(ora_lfc)

    # --- 3. Compute ORA for Each Comparison ---
    # Iterate through each dataset contrast to perform dataset-specific over-representation analysis
    for (comp_name in names(lfc_list)) {
      gse_id <- comparison_gse_map[comp_name]

      pvals <- pval_list[[comp_name]]
      lfcs <- lfc_list[[comp_name]]
      ses <- se_list[[comp_name]]

      # Ensure alignment to rownames(cts) by gene name if named
      all_genes <- rownames(cts)
      if (!is.null(names(pvals))) pvals <- pvals[all_genes]
      if (!is.null(names(lfcs))) lfcs <- lfcs[all_genes]
      if (!is.null(names(ses))) ses <- ses[all_genes]

      valid_mask <- !is.na(pvals) & !is.na(lfcs)
      universe_clusters <- all_genes[valid_mask]

      if (length(universe_clusters) == 0) next

      dataset_padj <- p.adjust(pvals[valid_mask], method = "BH")
      abs_lfc <- abs(lfcs[valid_mask])
      signed_lfc <- lfcs[valid_mask]

      universe_symbols <- unique(cluster_id_to_symbol[universe_clusters])
      universe_symbols <- universe_symbols[!is.na(universe_symbols) & universe_symbols != ""]

      if (length(universe_symbols) == 0) next

      # --- 3b. Build 2x2 tables for ALL globally valid terms ---
      # enricher() drops terms with zero DE genes, which would turn a = 0 into NA and bias the
      # meta-analysis upward. Count directly so zero-hit terms get a = 0 (strongly negative logOR).
      # Terms with no genes in this comparison's universe are untestable and remain NA.
      t2g_u <- term2gene[term2gene$gene %in% universe_symbols, ]
      bg_num <- as.numeric(table(factor(t2g_u$term, levels = valid_terms)))
      bg_den <- length(universe_symbols)
      testable <- bg_num > 0

      if (!any(testable)) next

      bg_num_sub <- bg_num[testable]
      valid_terms_sub <- valid_terms[testable]

      de_counts_log <- list()

      for (d in directions) {
        de_mask <- switch(d,
          "all"  = dataset_padj < ora_padj & abs_lfc >= abs_ora_lfc,
          "up"   = dataset_padj < ora_padj & signed_lfc >= abs_ora_lfc,
          "down" = dataset_padj < ora_padj & signed_lfc <= -abs_ora_lfc
        )
        de_clusters <- universe_clusters[de_mask]
        de_symbols <- unique(cluster_id_to_symbol[de_clusters])
        de_symbols <- de_symbols[!is.na(de_symbols) & de_symbols != ""]
        de_size <- length(de_symbols)
        de_counts_log[[d]] <- de_size

        # A comparison with 0 DE genes cannot be tested for enrichment and would yield
        # artifactual positive logOR values if included with Haldane-Anscombe correction (+0.5).
        if (de_size == 0) next

        # Record overlapping DE genes contributing to GO terms for this comparison
        overlap_mask <- t2g_u$gene %in% de_symbols
        if (any(overlap_mask)) {
          comp_overlaps <- t2g_u[overlap_mask, c("term", "gene")]
          comp_overlaps$comparison <- comp_name
          comp_overlaps$gse <- gse_id
          ora_contributing_genes[[d]][[comp_name]] <- comp_overlaps
        }

        a <- as.numeric(table(factor(t2g_u$term[overlap_mask], levels = valid_terms)))
        a_sub <- a[testable]
        b <- de_size - a_sub
        c_val <- bg_num_sub - a_sub
        d_val <- bg_den - de_size - c_val

        # Haldane-Anscombe correction (+0.5)
        a_adj <- a_sub + 0.5
        b_adj <- b + 0.5
        c_adj <- c_val + 0.5
        d_adj <- d_val + 0.5

        logOR <- log((a_adj * d_adj) / (b_adj * c_adj))
        se_logOR <- sqrt(1 / a_adj + 1 / b_adj + 1 / c_adj + 1 / d_adj)

        ora_logor_list[[d]][[comp_name]] <- setNames(logOR, valid_terms_sub)
        ora_se_list[[d]][[comp_name]] <- setNames(se_logOR, valid_terms_sub)
        ora_comparisons[[d]] <- c(ora_comparisons[[d]], comp_name)
        ora_gse_map[[d]][comp_name] <- gse_id
      }

      cat(sprintf(
        "  [%s] Eligible genes: %d | Mapped universe genes: %d | Mapped DE genes (all: %s, up: %s, down: %s)\n",
        comp_name, length(universe_clusters), length(universe_symbols),
        if (!is.null(de_counts_log[["all"]])) as.character(de_counts_log[["all"]]) else "-",
        if (!is.null(de_counts_log[["up"]])) as.character(de_counts_log[["up"]]) else "-",
        if (!is.null(de_counts_log[["down"]])) as.character(de_counts_log[["down"]]) else "-"
      ))
    }

    # --- 5. Aggregate Effect Sizes and Meta-Analyze Each Direction ---
    results_by_dir <- list()
    contrib_by_dir <- list()
    summary_by_dir <- list()
    debug_printed <- FALSE

    for (d in directions) {
      cat(sprintf("\n--- Meta-Analyzing GO Terms for Direction: %s ---\n", toupper(d)))

      curr_logor_list <- ora_logor_list[[d]]
      curr_se_list <- ora_se_list[[d]]
      curr_comparisons <- ora_comparisons[[d]]
      curr_gse_map <- ora_gse_map[[d]]

      if (length(curr_logor_list) == 0) {
        cat(sprintf("  - No valid ORA results computed for direction '%s'.\n", d))
        next
      }

      # Pre-aggregate all contributing DE gene overlaps for direction d
      all_overlaps_d <- if (length(ora_contributing_genes[[d]]) > 0) {
        do.call(rbind, ora_contributing_genes[[d]])
      } else {
        NULL
      }
      overlaps_by_term <- if (!is.null(all_overlaps_d) && nrow(all_overlaps_d) > 0) {
        split(all_overlaps_d, all_overlaps_d$term)
      } else {
        list()
      }

      all_go_ids <- unique(unlist(lapply(curr_logor_list, names)))

      ora_yi <- matrix(NA, nrow = length(all_go_ids), ncol = length(curr_logor_list), dimnames = list(all_go_ids, curr_comparisons))
      ora_vi <- matrix(NA, nrow = length(all_go_ids), ncol = length(curr_se_list), dimnames = list(all_go_ids, curr_comparisons))

      for (comp_name in curr_comparisons) {
        ora_yi[names(curr_logor_list[[comp_name]]), comp_name] <- curr_logor_list[[comp_name]]
        ora_vi[names(curr_se_list[[comp_name]]), comp_name] <- curr_se_list[[comp_name]]^2
      }

      cat(sprintf("  Meta-analyzing %d GO terms across %d comparisons (direction: %s)...\n",
                  length(all_go_ids), length(curr_comparisons), d))

      ora_meta_res <- list()
      for (i in seq_along(all_go_ids)) {
        go_id <- all_go_ids[i]
        yi_go <- ora_yi[i, ]
        vi_go <- ora_vi[i, ]

        # --- 6. Meta-Analyze Each GO Term ---
        # Fit 3-Level Random Effects (w/ clubSandwich RVE) and Fixed Effects models
        valid_go <- !is.na(yi_go) & !is.na(vi_go)
        valid_gse_go <- curr_gse_map[valid_go]
        valid_comp_go <- curr_comparisons[valid_go]
        n_valid_comps <- sum(valid_go)

        if (length(unique(valid_gse_go)) >= max(2, min_ds)) {
          # Extract contributing genes for this GO term across valid comparisons
          term_sub <- overlaps_by_term[[go_id]]
          if (!is.null(term_sub) && nrow(term_sub) > 0) {
            term_sub <- term_sub[term_sub$comparison %in% valid_comp_go, , drop = FALSE]
          }

          if (!is.null(term_sub) && nrow(term_sub) > 0) {
            gene_comp_counts <- tapply(term_sub$comparison, term_sub$gene, function(x) length(unique(x)))
            # Sort by frequency descending, then alphabetically
            sorted_genes <- names(gene_comp_counts)[order(-gene_comp_counts, names(gene_comp_counts))]
            sorted_counts <- as.integer(gene_comp_counts[sorted_genes])

            term_gene_count <- length(sorted_genes)
            term_genes <- paste(sorted_genes, collapse = "/")
            gene_freq_strs <- sprintf("%s (%d/%d)", sorted_genes, sorted_counts, n_valid_comps)
            term_genes_freq <- paste(gene_freq_strs, collapse = "; ")
          } else {
            term_gene_count <- 0L
            term_genes <- ""
            term_genes_freq <- ""
          }

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

          if (!debug_printed) {
            cat(sprintf("\n--- DEBUG: Covariance Matrices for GO term %s (direction: %s) ---\n", go_id, d))
            cat("Comparisons:\n")
            print(valid_comp_go)
            cat("\nClusters (valid_cov_cluster):\n")
            print(valid_cov_cluster)
            cat("\nCorrelation Matrix (derived from V_mat_logor):\n")
            print(cov2cor(V_mat_logor))
            cat("-----------------------------------------------\n\n")
            debug_printed <- TRUE
          }

          # LogOR Meta-Analysis
          res_re <- tryCatch(
            metafor::rma.mv(yi = yi_go[valid_go], V = V_mat_logor, random = ~ 1 | valid_gse_go / valid_comp_go, method = "REML"),
            error = function(e) {
              cat(sprintf("  - [logOR %s] 3-level RE model failed for %s: %s\n", d, go_id, e$message))
              NULL
            }
          )

          if (!is.null(res_re)) {
            res_re_rob <- tryCatch(
              metafor::robust(res_re, cluster = valid_gse_go, clubSandwich = TRUE),
              error = function(e) {
                cat(sprintf("  - [logOR %s] clubSandwich RVE failed for %s (using model SE): %s\n", d, go_id, e$message))
                NULL
              }
            )
            if (!is.null(res_re_rob)) {
              res_re <- res_re_rob
            }
          }

          res_fe <- tryCatch(
            metafor::rma.mv(yi = yi_go[valid_go], V = V_mat_logor),
            error = function(e) {
              cat(sprintf("  - [logOR FE %s] model failed for %s: %s\n", d, go_id, e$message))
              NULL
            }
          )

          if (!is.null(res_re) || !is.null(res_fe)) {
            re_est <- if (!is.null(res_re)) as.numeric(res_re$beta)[1] else NA_real_
            re_se  <- if (!is.null(res_re)) as.numeric(res_re$se)[1] else NA_real_
            re_p   <- if (!is.null(res_re)) one_sided_p(re_est, as.numeric(res_re$pval)[1]) else NA_real_

            fe_est <- if (!is.null(res_fe)) as.numeric(res_fe$beta)[1] else NA_real_
            fe_se  <- if (!is.null(res_fe)) as.numeric(res_fe$se)[1] else NA_real_
            fe_p   <- if (!is.null(res_fe)) one_sided_p(fe_est, as.numeric(res_fe$pval)[1]) else NA_real_

            ora_meta_res[[go_id]] <- data.frame(
              ID = go_id,
              Ontology = term2ontology$ontology[match(go_id, term2ontology$term)],
              Description = term2name$name[match(go_id, term2name$term)],
              Direction = d,
              logOR_RE_estimate = re_est,
              logOR_RE_SE = re_se,
              logOR_RE_pval = re_p,
              logOR_FE_estimate = fe_est,
              logOR_FE_SE = fe_se,
              logOR_FE_pval = fe_p,
              N_comparisons = n_valid_comps,
              N_datasets = length(unique(valid_gse_go)),
              Gene_Count = term_gene_count,
              Genes = term_genes,
              Genes_Frequency = term_genes_freq,
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

        # Reorder columns to put padj next to pval, followed by gene info
        ora_meta_df <- ora_meta_df[, c(
          "ID", "Ontology", "Description", "Direction",
          "logOR_RE_estimate", "logOR_RE_SE", "logOR_RE_pval", "logOR_RE_padj",
          "logOR_FE_estimate", "logOR_FE_SE", "logOR_FE_pval", "logOR_FE_padj",
          "N_comparisons", "N_datasets",
          "Gene_Count", "Genes", "Genes_Frequency"
        )]

        ora_meta_df <- ora_meta_df[order(ora_meta_df$logOR_RE_pval, na.last = TRUE), ]

        # Save direction-specific CSV
        out_file_d <- file.path(out_dir, sprintf("ORA_Meta_Analysis_Results_%s.csv", d))
        write.csv(ora_meta_df, out_file_d, row.names = FALSE)
        cat(sprintf("  - Saved ORA Meta-Analysis Results (%s) to %s\n", d, out_file_d))

        # For "all", also write to standard ORA_Meta_Analysis_Results.csv for backward compatibility
        if (d == "all") {
          out_file_legacy <- file.path(out_dir, "ORA_Meta_Analysis_Results.csv")
          write.csv(ora_meta_df, out_file_legacy, row.names = FALSE)
          cat(sprintf("  - Saved ORA Meta-Analysis Results (legacy default) to %s\n", out_file_legacy))
        }

        results_by_dir[[d]] <- ora_meta_df

        # --- 7b. Export Detailed Gene Contributions and Gene Pleiotropy Summaries ---
        if (!is.null(all_overlaps_d) && nrow(all_overlaps_d) > 0) {
          # Keep only overlaps for terms evaluated in the meta-analysis
          overlaps_tested <- all_overlaps_d[all_overlaps_d$term %in% ora_meta_df$ID, , drop = FALSE]

          if (nrow(overlaps_tested) > 0) {
            term_valid_n <- setNames(ora_meta_df$N_comparisons, ora_meta_df$ID)
            tg_key <- paste(overlaps_tested$term, overlaps_tested$gene, sep = "@@@")
            tg_splits <- split(overlaps_tested, tg_key)

            contrib_rows <- lapply(tg_splits, function(df_tg) {
              t_id <- df_tg$term[1]
              g_id <- df_tg$gene[1]
              comps <- unique(df_tg$comparison)
              gses  <- unique(df_tg$gse)
              n_comps <- length(comps)
              tot_comps <- if (!is.na(term_valid_n[t_id])) term_valid_n[t_id] else length(curr_comparisons)
              pct_comps <- round(100 * n_comps / tot_comps, 1)

              data.frame(
                GO_ID = t_id,
                Ontology = term2ontology$ontology[match(t_id, term2ontology$term)],
                Description = term2name$name[match(t_id, term2name$term)],
                Direction = d,
                Gene = g_id,
                N_Comparisons = n_comps,
                Total_Comparisons = tot_comps,
                Percent_Comparisons = pct_comps,
                N_Datasets = length(gses),
                Comparisons = paste(sort(comps), collapse = "; "),
                Datasets = paste(sort(gses), collapse = "; "),
                stringsAsFactors = FALSE
              )
            })
            contrib_df <- do.call(rbind, contrib_rows)
            contrib_df <- contrib_df[order(contrib_df$GO_ID, -contrib_df$N_Comparisons, contrib_df$Gene), ]

            contrib_file_d <- file.path(out_dir, sprintf("ORA_Meta_Analysis_Gene_Contributions_%s.csv", d))
            write.csv(contrib_df, contrib_file_d, row.names = FALSE)
            cat(sprintf("  - Saved Detailed Gene Contributions (%s) to %s\n", d, contrib_file_d))

            if (d == "all") {
              contrib_file_legacy <- file.path(out_dir, "ORA_Meta_Analysis_Gene_Contributions.csv")
              write.csv(contrib_df, contrib_file_legacy, row.names = FALSE)
              cat(sprintf("  - Saved Detailed Gene Contributions (legacy default) to %s\n", contrib_file_legacy))
            }
            contrib_by_dir[[d]] <- contrib_df

            # Export gene-level summary across significant terms (padj < report_padj)
            sig_term_ids <- ora_meta_df$ID[!is.na(ora_meta_df$logOR_RE_padj) & ora_meta_df$logOR_RE_padj < report_padj]
            sig_overlaps <- overlaps_tested[overlaps_tested$term %in% sig_term_ids, , drop = FALSE]

            if (nrow(sig_overlaps) > 0) {
              g_splits <- split(sig_overlaps, sig_overlaps$gene)
              gene_summary_rows <- lapply(g_splits, function(df_g) {
                g_id <- df_g$gene[1]
                terms <- unique(df_g$term)
                term_names <- term2name$name[match(terms, term2name$term)]
                term_desc <- sprintf("%s (%s)", terms, term_names)
                comps <- unique(df_g$comparison)
                gses <- unique(df_g$gse)
                data.frame(
                  Gene = g_id,
                  Direction = d,
                  N_Significant_Terms = length(terms),
                  N_Comparisons = length(comps),
                  N_Datasets = length(gses),
                  Significant_Terms = paste(term_desc, collapse = "; "),
                  Comparisons = paste(sort(comps), collapse = "; "),
                  Datasets = paste(sort(gses), collapse = "; "),
                  stringsAsFactors = FALSE
                )
              })
              gene_summary_df <- do.call(rbind, gene_summary_rows)
              gene_summary_df <- gene_summary_df[order(-gene_summary_df$N_Significant_Terms, -gene_summary_df$N_Comparisons, gene_summary_df$Gene), ]

              gene_sum_file_d <- file.path(out_dir, sprintf("ORA_Meta_Analysis_Gene_Summary_%s.csv", d))
              write.csv(gene_summary_df, gene_sum_file_d, row.names = FALSE)
              cat(sprintf("  - Saved Gene-Level Summary across significant terms (%s) to %s\n", d, gene_sum_file_d))

              if (d == "all") {
                gene_sum_file_legacy <- file.path(out_dir, "ORA_Meta_Analysis_Gene_Summary.csv")
                write.csv(gene_summary_df, gene_sum_file_legacy, row.names = FALSE)
                cat(sprintf("  - Saved Gene-Level Summary (legacy default) to %s\n", gene_sum_file_legacy))
              }
              summary_by_dir[[d]] <- gene_summary_df
            }
          }
        }

        # --- 8. Reporting-only semantic clustering (after BH; full table above is unchanged) ---
        clustered <- tryCatch(
          cluster_significant_go_terms(ora_meta_df, sig_padj = report_padj, sim_cutoff = sem_sim_cutoff),
          error = function(e) {
            cat(sprintf("  - [%s] GOSemSim clustering failed: %s\n", d, e$message))
            NULL
          }
        )
        if (!is.null(clustered)) {
          clustered <- clustered[order(clustered$Ontology, clustered$Cluster, !clustered$Is_Representative, clustered$logOR_RE_padj), ]
          clusters_file <- file.path(out_dir, sprintf("ORA_Meta_Analysis_Significant_Term_Clusters_%s.csv", d))
          write.csv(clustered, clusters_file, row.names = FALSE)
          reps <- clustered[clustered$Is_Representative, ]
          reps <- reps[order(reps$logOR_RE_padj), ]
          reps_file <- file.path(out_dir, sprintf("ORA_Meta_Analysis_Representative_Terms_%s.csv", d))
          write.csv(reps, reps_file, row.names = FALSE)
          cat(sprintf("  - %d significant terms (padj < %g) -> %d clusters (Wang sim >= %g). Saved clusters to %s and representatives to %s\n",
                      nrow(clustered), report_padj, nrow(reps), sem_sim_cutoff, clusters_file, reps_file))
        } else {
          cat(sprintf("  - [%s] No significant terms (padj < %g) to cluster.\n", d, report_padj))
        }
      } else {
        cat(sprintf("  - [%s] Not enough valid datasets/comparisons to run meta-analysis on any GO term.\n", d))
      }
    }

    # If multiple directions were analyzed, also save a combined CSV table
    if (length(results_by_dir) > 1) {
      combined_df <- do.call(rbind, results_by_dir)
      combined_file <- file.path(out_dir, "ORA_Meta_Analysis_Results_combined.csv")
      write.csv(combined_df, combined_file, row.names = FALSE)
      cat(sprintf("  - Saved Combined ORA Meta-Analysis Results across all directions to %s\n", combined_file))

      if (length(contrib_by_dir) > 0) {
        combined_contrib <- do.call(rbind, contrib_by_dir)
        combined_contrib_file <- file.path(out_dir, "ORA_Meta_Analysis_Gene_Contributions_combined.csv")
        write.csv(combined_contrib, combined_contrib_file, row.names = FALSE)
        cat(sprintf("  - Saved Combined Detailed Gene Contributions across all directions to %s\n", combined_contrib_file))
      }

      if (length(summary_by_dir) > 0) {
        combined_summary <- do.call(rbind, summary_by_dir)
        combined_sum_file <- file.path(out_dir, "ORA_Meta_Analysis_Gene_Summary_combined.csv")
        write.csv(combined_summary, combined_sum_file, row.names = FALSE)
        cat(sprintf("  - Saved Combined Gene-Level Summary across all directions to %s\n", combined_sum_file))
      }
    }

    if (length(results_by_dir) == 0) {
      cat("  - No ORA meta-analysis results could be computed for any direction.\n")
    }

    return(invisible(results_by_dir))
  } else {
    cat("  - GAF file not found or not provided. Skipping ORA meta-analysis.\n")
  }
}
