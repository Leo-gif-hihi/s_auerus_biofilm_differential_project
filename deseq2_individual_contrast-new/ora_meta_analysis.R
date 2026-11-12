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

# Ensure required Bioconductor and CRAN packages are installed without polluting the search path.
.ensure_ora_packages <- function(bioc_pkgs = character(0), cran_pkgs = character(0)) {
  for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
      requireNamespace(pkg, quietly = TRUE)
    }
  }
  for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      requireNamespace(pkg, quietly = TRUE)
    }
  }
}

# One-sided p-value testing estimate > 0 from a two-sided p-value.
# Positive estimate: halve p. Negative estimate: 1 - p/2 (never significant).
.one_sided_p <- function(est, p) ifelse(est > 0, p / 2, 1 - p / 2)

# Map each cluster_id to the first token in All_Genes_Locus_Tags matching target_genes.
.map_clusters_by_tokens <- function(cluster_annotation, target_genes) {
  mapped <- sapply(cluster_annotation$All_Genes_Locus_Tags, function(tags) {
    if (is.na(tags) || tags == "") return(NA_character_)
    toks <- strsplit(as.character(tags), "/")[[1]]
    matches <- toks[toks %in% target_genes]
    if (length(matches) > 0) matches[1] else NA_character_
  })
  names(mapped) <- cluster_annotation$cluster_id
  mapped
}

# Apply global gene-set size filtering against the study universe so tested terms
# are uniform across comparisons.
.filter_term2gene_universe <- function(term2gene, cluster_id_to_gene, study_clusters, min_gs_size, max_gs_size) {
  global_genes <- unique(cluster_id_to_gene[study_clusters])
  global_genes <- global_genes[!is.na(global_genes) & global_genes != ""]
  if (length(global_genes) == 0) {
    return(list(term2gene = term2gene[0, , drop = FALSE], valid_terms = character(0), global_genes = character(0)))
  }
  t2g_study <- term2gene[term2gene$gene %in% global_genes, , drop = FALSE]
  term_sizes <- table(t2g_study$term)
  valid_terms <- names(term_sizes)[term_sizes >= min_gs_size & term_sizes <= max_gs_size]
  t2g_study <- t2g_study[t2g_study$term %in% valid_terms, , drop = FALSE]
  list(term2gene = t2g_study, valid_terms = valid_terms, global_genes = global_genes)
}

# Parse GAF annotations, propagate GO DAG ancestors, and filter to study universe.
.prepare_go_annotations <- function(gaf_file, cluster_annotation, study_clusters, exclude_evidence = c("ND"), min_gs_size = 5, max_gs_size = 500) {
  .ensure_ora_packages(bioc_pkgs = c("clusterProfiler", "GO.db"))

  gaf <- read.delim(gaf_file, comment.char = "!", header = FALSE, stringsAsFactors = FALSE, quote = "")
  n_raw <- nrow(gaf)
  # Exclude negative qualifiers in V4 and uninformative evidence codes in V7
  gaf <- gaf[!grepl("(^|\\|)NOT(\\||$)", gaf$V4), ]
  if (!is.null(exclude_evidence) && length(exclude_evidence) > 0) {
    gaf <- gaf[!gaf$V7 %in% exclude_evidence, ]
  }
  cat(sprintf("  GAF filtering: %d raw -> %d annotations retained (filtered NOT qualifiers; excluded evidence: %s)\n",
              n_raw, nrow(gaf), if (length(exclude_evidence) > 0) paste(exclude_evidence, collapse = ", ") else "none"))

  term2gene <- unique(data.frame(term = gaf$V5, gene = gaf$V3, stringsAsFactors = FALSE))
  n_direct <- nrow(term2gene)
  term2gene <- clusterProfiler::buildGOmap(term2gene)
  term2gene <- unique(data.frame(term = as.character(term2gene[[1]]), gene = as.character(term2gene[[2]]), stringsAsFactors = FALSE))
  # Drop the ontology root terms (BP/MF/CC), which carry no specific signal
  term2gene <- term2gene[grepl("^GO:", term2gene$term) & !term2gene$term %in% c("GO:0008150", "GO:0003674", "GO:0005575"), ]
  cat(sprintf("  Propagated GO annotations: %d direct -> %d direct+indirect term-gene pairs\n", n_direct, nrow(term2gene)))

  annotated_genes <- unique(gaf$V3)
  cluster_id_to_gene <- .map_clusters_by_tokens(cluster_annotation, annotated_genes)

  filt <- .filter_term2gene_universe(term2gene, cluster_id_to_gene, study_clusters, min_gs_size, max_gs_size)
  cat(sprintf("  Global gene-set size filtering (%d-%d): %d GO terms retained across %d study genes\n",
              min_gs_size, max_gs_size, length(filt$valid_terms), length(filt$global_genes)))

  goterms <- as.data.frame(GO.db::GOTERM)
  term2name <- unique(data.frame(term = goterms$go_id, name = goterms$Term, stringsAsFactors = FALSE))
  term2ontology <- unique(data.frame(term = goterms$go_id, ontology = goterms$Ontology, stringsAsFactors = FALSE))

  list(
    term2gene = filt$term2gene,
    term2name = term2name,
    term2ontology = term2ontology,
    valid_terms = filt$valid_terms,
    global_genes = filt$global_genes,
    cluster_id_to_gene = cluster_id_to_gene,
    gene_to_symbol = NULL
  )
}

# Fetch and cache KEGG pathway annotations for the target organism (default: NCTC 8325 = "sao").
.load_kegg_annotations <- function(organism = "sao", cache_file = "kegg_sao_cache.rds", refresh = FALSE) {
  .ensure_ora_packages(bioc_pkgs = c("clusterProfiler"))

  if (!isTRUE(refresh) && !is.null(cache_file) && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(cached) && is.list(cached) &&
        all(c("term2gene", "term2name") %in% names(cached)) &&
        identical(cached$organism, organism)) {
      cat(sprintf("  Loaded cached KEGG annotations (%s) from %s\n", organism, cache_file))
      return(cached)
    }
  }

  cat(sprintf("  Downloading KEGG annotations for organism '%s'...\n", organism))
  kegg <- clusterProfiler::download_KEGG(organism, keggType = "KEGG", keyType = "kegg")
  term2gene <- unique(data.frame(
    term = as.character(kegg$KEGGPATHID2EXTID$from),
    gene = as.character(kegg$KEGGPATHID2EXTID$to),
    stringsAsFactors = FALSE
  ))
  term2name <- unique(data.frame(
    term = as.character(kegg$KEGGPATHID2NAME$from),
    name = as.character(kegg$KEGGPATHID2NAME$to),
    stringsAsFactors = FALSE
  ))
  # Strip organism suffix appended by KEGG (e.g., " - Staphylococcus aureus subsp. aureus NCTC 8325")
  term2name$name <- sub(" - Staphylococcus aureus.*$", "", term2name$name)

  res <- list(organism = organism, term2gene = term2gene, term2name = term2name)
  if (!is.null(cache_file) && nzchar(cache_file)) {
    cache_dir <- dirname(cache_file)
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    tryCatch(saveRDS(res, cache_file), error = function(e) {
      cat("  ⚠ Could not save KEGG cache to", cache_file, ":", e$message, "\n")
    })
    cat(sprintf("  Cached KEGG annotations (%s) to %s\n", organism, cache_file))
  }
  res
}

# Build hybrid cluster_id -> KEGG locus_tag (SAOUHSC_*) mapping and readable gene symbols.
.prepare_kegg_annotations <- function(cluster_annotation, study_clusters, organism = "sao", kegg_cache_file = "kegg_sao_cache.rds", refresh_kegg = FALSE, operon_mapping_file = "operon_gene_mapping.csv", min_gs_size = 5, max_gs_size = 500, exclude_overview_011 = TRUE) {
  kegg_raw <- .load_kegg_annotations(organism = organism, cache_file = kegg_cache_file, refresh = refresh_kegg)
  term2gene <- kegg_raw$term2gene
  term2name <- kegg_raw$term2name

  if (isTRUE(exclude_overview_011)) {
    overview_pat <- paste0("^", organism, "011[0-9]{2}$")
    n_before <- length(unique(term2gene$term))
    term2gene <- term2gene[!grepl(overview_pat, term2gene$term), , drop = FALSE]
    term2name <- term2name[!grepl(overview_pat, term2name$term), , drop = FALSE]
    cat(sprintf("  KEGG filtering: excluded 011xx global overview maps (%d -> %d pathways)\n",
                n_before, length(unique(term2gene$term))))
  }

  kegg_genes <- unique(term2gene$gene)

  # 1. Direct token match against All_Genes_Locus_Tags
  cluster_id_to_gene <- .map_clusters_by_tokens(cluster_annotation, kegg_genes)
  n_direct <- sum(!is.na(cluster_id_to_gene[study_clusters]))

  # 2. Fallback via operon_gene_mapping.csv (final_cluster_id -> locus_tag) for ortholog-mapped clusters
  op_df <- NULL
  if (!is.null(operon_mapping_file) && file.exists(operon_mapping_file)) {
    op_df <- tryCatch(read.csv(operon_mapping_file, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(op_df) && all(c("final_cluster_id", "locus_tag") %in% colnames(op_df))) {
      op_sub <- op_df[!is.na(op_df$final_cluster_id) & op_df$final_cluster_id != "" &
                        !is.na(op_df$locus_tag) & op_df$locus_tag %in% kegg_genes, , drop = FALSE]
      op_sub <- op_sub[!duplicated(op_sub$final_cluster_id), , drop = FALSE]
      missing_cids <- names(cluster_id_to_gene)[is.na(cluster_id_to_gene)]
      hit <- op_sub[op_sub$final_cluster_id %in% missing_cids, , drop = FALSE]
      if (nrow(hit) > 0) {
        cluster_id_to_gene[hit$final_cluster_id] <- hit$locus_tag
      }
    }
  }
  n_total <- sum(!is.na(cluster_id_to_gene[study_clusters]))
  cat(sprintf("  KEGG gene mapping (%s): %d direct + %d operon-ortholog fallback = %d / %d study clusters mapped\n",
              organism, n_direct, n_total - n_direct, n_total, length(study_clusters)))

  # Build readable gene symbol lookup for SAOUHSC_* locus tags (operon_gene_mapping + Representative_Gene)
  gene_to_symbol <- setNames(rep(NA_character_, length(kegg_genes)), kegg_genes)
  if (!is.null(op_df) && all(c("locus_tag", "gene_symbol") %in% colnames(op_df))) {
    sym_rows <- op_df[!is.na(op_df$locus_tag) & op_df$locus_tag %in% kegg_genes &
                        !is.na(op_df$gene_symbol) & op_df$gene_symbol != "" &
                        op_df$gene_symbol != "-" & op_df$gene_symbol != op_df$locus_tag, , drop = FALSE]
    sym_rows <- sym_rows[!duplicated(sym_rows$locus_tag), , drop = FALSE]
    if (nrow(sym_rows) > 0) {
      gene_to_symbol[sym_rows$locus_tag] <- sym_rows$gene_symbol
    }
  }
  if ("Representative_Gene" %in% colnames(cluster_annotation)) {
    mapped_idx <- which(!is.na(cluster_id_to_gene))
    for (idx in mapped_idx) {
      kg <- cluster_id_to_gene[idx]
      rg <- as.character(cluster_annotation$Representative_Gene[idx])
      if (is.na(gene_to_symbol[kg]) && !is.na(rg) && rg != "" && rg != "-" && rg != kg && !grepl("^SAOUHSC_", rg)) {
        gene_to_symbol[kg] <- rg
      }
    }
  }

  filt <- .filter_term2gene_universe(term2gene, cluster_id_to_gene, study_clusters, min_gs_size, max_gs_size)
  cat(sprintf("  Global gene-set size filtering (%d-%d): %d KEGG pathways retained across %d study genes\n",
              min_gs_size, max_gs_size, length(filt$valid_terms), length(filt$global_genes)))

  term2ontology <- data.frame(term = filt$valid_terms, ontology = "KEGG", stringsAsFactors = FALSE)

  list(
    term2gene = filt$term2gene,
    term2name = term2name,
    term2ontology = term2ontology,
    valid_terms = filt$valid_terms,
    global_genes = filt$global_genes,
    cluster_id_to_gene = cluster_id_to_gene,
    gene_to_symbol = gene_to_symbol
  )
}

# Build measured-universe operon structures (Option B) shared by Meta-ORA and FE-ORA.
.setup_operon_units <- function(operon_mapping_file, study_clusters, cluster_id_to_gene, term2gene, operon_de_threshold = NULL) {
  empty <- list(enabled = FALSE, op_to_cids = list(), cid_to_ops = list(), op_labels = character(0), term_to_ops = list())
  if (is.null(operon_mapping_file) || !file.exists(operon_mapping_file)) return(empty)

  cat("Loading Operon mapping file:", operon_mapping_file, "\n")
  op_map <- tryCatch(
    read.csv(operon_mapping_file, stringsAsFactors = FALSE),
    error = function(e) {
      cat("  ⚠ Failed to read operon mapping file:", e$message, "\n")
      NULL
    }
  )
  if (is.null(op_map) || !all(c("final_cluster_id", "operon_id") %in% colnames(op_map))) return(empty)

  valid_op <- op_map[!is.na(op_map$final_cluster_id) & op_map$final_cluster_id != "", , drop = FALSE]
  if (nrow(valid_op) == 0) return(empty)

  op_to_cids <- lapply(split(valid_op$final_cluster_id, valid_op$operon_id), unique)
  cid_to_ops <- lapply(split(valid_op$operon_id, valid_op$final_cluster_id), unique)

  op_labels <- character(0)
  for (oid in unique(valid_op$operon_id)) {
    sub_op <- valid_op[valid_op$operon_id == oid, , drop = FALSE]
    syms <- sub_op$gene_symbol[!is.na(sub_op$gene_symbol) & sub_op$gene_symbol != ""]
    sym_str <- if (length(syms) > 0) paste(unique(syms), collapse = "/") else ""
    op_labels[as.character(oid)] <- if (sym_str != "") paste0("Operon_", oid, " (", sym_str, ")") else paste0("Operon_", oid)
  }

  unassigned <- setdiff(study_clusters, names(cid_to_ops))
  if (length(unassigned) > 0) {
    for (u_cid in unassigned) {
      fake_id <- paste0("Singleton_", u_cid)
      op_to_cids[[fake_id]] <- u_cid
      cid_to_ops[[u_cid]] <- fake_id
      op_labels[fake_id] <- fake_id
    }
  }

  gene_to_terms <- split(term2gene$term, term2gene$gene)
  cid_to_terms <- lapply(study_clusters, function(cid) {
    sym <- cluster_id_to_gene[[cid]]
    if (!is.na(sym) && sym != "" && sym %in% names(gene_to_terms)) {
      unique(gene_to_terms[[sym]])
    } else {
      character(0)
    }
  })
  names(cid_to_terms) <- study_clusters

  op_to_terms <- lapply(op_to_cids, function(cids) {
    unique(unlist(cid_to_terms[intersect(cids, study_clusters)]))
  })

  op_names_vec <- rep(names(op_to_terms), lengths(op_to_terms))
  term_vals_vec <- unlist(op_to_terms, use.names = FALSE)
  term_to_ops <- if (length(term_vals_vec) > 0) {
    lapply(split(op_names_vec, term_vals_vec), unique)
  } else {
    list()
  }

  if (!is.null(operon_de_threshold)) {
    cat(sprintf("  ✓ Operon dual-ORA enabled: %d operon units across %d study clusters (Option B, threshold >= %.2f).\n",
                length(op_to_cids), length(study_clusters), operon_de_threshold))
  }

  list(enabled = TRUE, op_to_cids = op_to_cids, cid_to_ops = cid_to_ops, op_labels = op_labels, term_to_ops = term_to_ops)
}

# Format gene symbols corresponding to a character vector of canonical gene IDs.
.resolve_gene_symbols <- function(genes, gene_to_symbol) {
  if (is.null(gene_to_symbol) || length(genes) == 0) return(genes)
  syms <- unname(gene_to_symbol[genes])
  ifelse(!is.na(syms) & syms != "", syms, genes)
}

# Format gene label for frequency strings (e.g. "SAOUHSC_01807 [pyk]" when symbol is available).
.format_gene_freq_labels <- function(genes, gene_to_symbol) {
  if (is.null(gene_to_symbol) || length(genes) == 0) return(genes)
  syms <- unname(gene_to_symbol[genes])
  ifelse(!is.na(syms) & syms != "", sprintf("%s [%s]", genes, syms), genes)
}

# Shared per-comparison 2x2 Haldane-Anscombe ORA + 3-level RE / FE meta-analysis engine.
.run_geneset_meta_ora_core <- function(prep, lfc_list, pval_list, se_list, cts, comparison_gse_map, comparison_control_map = NULL, min_ds, out_dir, ora_padj = 0.05, ora_lfc = 1.0, directions = c("all", "up", "down"), sem_sim_cutoff = 0.7, report_padj = 0.05, operon_mapping_file = "operon_gene_mapping.csv", operon_de_threshold = 0.5, file_prefix = "ORA_Meta_Analysis", id_col_name = "GO_ID", term_noun = "GO terms", write_legacy_all = TRUE, run_clustering = TRUE) {
  .ensure_ora_packages(cran_pkgs = c("clubSandwich"))

  term2gene <- prep$term2gene
  term2name <- prep$term2name
  term2ontology <- prep$term2ontology
  valid_terms <- prep$valid_terms
  cluster_id_to_symbol <- prep$cluster_id_to_gene
  gene_to_symbol <- prep$gene_to_symbol
  has_symbol_map <- !is.null(gene_to_symbol)

  if (length(valid_terms) == 0) {
    cat(sprintf("  - No valid %s retained after size filtering.\n", term_noun))
    return(invisible(list()))
  }

  all_genes <- rownames(cts)
  op_info <- .setup_operon_units(operon_mapping_file, all_genes, cluster_id_to_symbol, term2gene, operon_de_threshold = operon_de_threshold)
  operon_enabled <- op_info$enabled
  op_to_cids <- op_info$op_to_cids
  op_labels <- op_info$op_labels
  term_to_ops <- op_info$term_to_ops

  ora_logor_list <- setNames(lapply(directions, function(x) list()), directions)
  ora_se_list <- setNames(lapply(directions, function(x) list()), directions)
  ora_comparisons <- setNames(lapply(directions, function(x) character(0)), directions)
  ora_gse_map <- setNames(lapply(directions, function(x) character(0)), directions)
  ora_contributing_genes <- setNames(lapply(directions, function(x) list()), directions)

  ora_op_logor_list <- setNames(lapply(directions, function(x) list()), directions)
  ora_op_se_list <- setNames(lapply(directions, function(x) list()), directions)
  ora_contributing_operons <- setNames(lapply(directions, function(x) list()), directions)

  abs_ora_lfc <- abs(ora_lfc)

  # --- Compute ORA for Each Comparison ---
  for (comp_name in names(lfc_list)) {
    gse_id <- comparison_gse_map[comp_name]

    pvals <- pval_list[[comp_name]]
    lfcs <- lfc_list[[comp_name]]
    ses <- se_list[[comp_name]]

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

    # Build 2x2 tables for all globally valid terms so zero-hit terms are kept (a = 0)
    t2g_u <- term2gene[term2gene$gene %in% universe_symbols, , drop = FALSE]
    bg_num <- as.numeric(table(factor(t2g_u$term, levels = valid_terms)))
    bg_den <- length(universe_symbols)
    testable <- bg_num > 0
    if (!any(testable)) next

    bg_num_sub <- bg_num[testable]
    valid_terms_sub <- valid_terms[testable]

    de_counts_log <- list()
    de_op_counts_log <- list()

    if (operon_enabled) {
      op_meas_counts <- sapply(op_to_cids, function(cids) sum(cids %in% universe_clusters))
      valid_comp_ops <- names(op_to_cids)[op_meas_counts > 0]
      total_bg_ops <- length(valid_comp_ops)

      bg_op_num <- sapply(valid_terms_sub, function(t_id) {
        t_ops <- term_to_ops[[t_id]]
        if (is.null(t_ops)) 0L else sum(valid_comp_ops %in% t_ops)
      })
    }

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

      if (de_size == 0) next

      overlap_mask <- t2g_u$gene %in% de_symbols
      if (any(overlap_mask)) {
        comp_overlaps <- t2g_u[overlap_mask, c("term", "gene"), drop = FALSE]
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

      # --- Operon-Level ORA for this comparison (Option B: Measured Gene Universe) ---
      if (operon_enabled && total_bg_ops > 0) {
        op_de_counts <- sapply(valid_comp_ops, function(oid) {
          cids <- op_to_cids[[oid]]
          sum(cids %in% de_clusters)
        })

        op_is_de <- (op_de_counts / op_meas_counts[valid_comp_ops]) >= operon_de_threshold
        de_operons <- valid_comp_ops[op_is_de]
        de_op_size <- length(de_operons)
        de_op_counts_log[[d]] <- de_op_size

        if (de_op_size > 0) {
          a_op <- sapply(valid_terms_sub, function(t_id) {
            t_ops <- term_to_ops[[t_id]]
            if (is.null(t_ops)) 0L else sum(de_operons %in% t_ops)
          })
          b_op <- de_op_size - a_op
          c_op <- bg_op_num - a_op
          d_op <- total_bg_ops - de_op_size - c_op

          a_op_adj <- a_op + 0.5
          b_op_adj <- b_op + 0.5
          c_op_adj <- c_op + 0.5
          d_op_adj <- d_op + 0.5

          logor_op <- log((a_op_adj * d_op_adj) / (b_op_adj * c_op_adj))
          se_op <- sqrt(1 / a_op_adj + 1 / b_op_adj + 1 / c_op_adj + 1 / d_op_adj)

          ora_op_logor_list[[d]][[comp_name]] <- setNames(logor_op, valid_terms_sub)
          ora_op_se_list[[d]][[comp_name]] <- setNames(se_op, valid_terms_sub)

          hit_terms <- valid_terms_sub[a_op > 0]
          if (length(hit_terms) > 0) {
            comp_op_overlaps_list <- lapply(hit_terms, function(t_id) {
              t_ops <- term_to_ops[[t_id]]
              hit_ops <- intersect(de_operons, t_ops)
              if (length(hit_ops) > 0) {
                data.frame(
                  term = t_id,
                  operon_id = hit_ops,
                  operon_label = unname(op_labels[hit_ops]),
                  comparison = comp_name,
                  gse = gse_id,
                  stringsAsFactors = FALSE
                )
              } else {
                NULL
              }
            })
            comp_op_overlaps <- do.call(rbind, comp_op_overlaps_list)
            if (!is.null(comp_op_overlaps) && nrow(comp_op_overlaps) > 0) {
              ora_contributing_operons[[d]][[comp_name]] <- comp_op_overlaps
            }
          }
        }
      }
    }

    cat(sprintf(
      "  [%s] Eligible genes: %d | Mapped universe genes: %d | Mapped DE genes (all: %s, up: %s, down: %s)%s\n",
      comp_name, length(universe_clusters), length(universe_symbols),
      if (!is.null(de_counts_log[["all"]])) as.character(de_counts_log[["all"]]) else "-",
      if (!is.null(de_counts_log[["up"]])) as.character(de_counts_log[["up"]]) else "-",
      if (!is.null(de_counts_log[["down"]])) as.character(de_counts_log[["down"]]) else "-",
      if (operon_enabled) sprintf(" | Mapped DE operons (all: %s, up: %s, down: %s)",
        if (!is.null(de_op_counts_log[["all"]])) as.character(de_op_counts_log[["all"]]) else "-",
        if (!is.null(de_op_counts_log[["up"]])) as.character(de_op_counts_log[["up"]]) else "-",
        if (!is.null(de_op_counts_log[["down"]])) as.character(de_op_counts_log[["down"]]) else "-"
      ) else ""
    ))
  }

  # --- Aggregate Effect Sizes and Meta-Analyze Each Direction ---
  results_by_dir <- list()
  contrib_by_dir <- list()
  summary_by_dir <- list()
  debug_printed <- FALSE

  for (d in directions) {
    cat(sprintf("\n--- Meta-Analyzing %s for Direction: %s ---\n", term_noun, toupper(d)))

    curr_logor_list <- ora_logor_list[[d]]
    curr_se_list <- ora_se_list[[d]]
    curr_comparisons <- ora_comparisons[[d]]
    curr_gse_map <- ora_gse_map[[d]]

    if (length(curr_logor_list) == 0) {
      cat(sprintf("  - No valid ORA results computed for direction '%s'.\n", d))
      next
    }

    all_overlaps_d <- if (length(ora_contributing_genes[[d]]) > 0) do.call(rbind, ora_contributing_genes[[d]]) else NULL
    overlaps_by_term <- if (!is.null(all_overlaps_d) && nrow(all_overlaps_d) > 0) split(all_overlaps_d, all_overlaps_d$term) else list()

    all_op_overlaps_d <- if (operon_enabled && length(ora_contributing_operons[[d]]) > 0) do.call(rbind, ora_contributing_operons[[d]]) else NULL
    op_overlaps_by_term <- if (!is.null(all_op_overlaps_d) && nrow(all_op_overlaps_d) > 0) split(all_op_overlaps_d, all_op_overlaps_d$term) else list()

    all_term_ids <- unique(unlist(lapply(curr_logor_list, names)))

    ora_yi <- matrix(NA, nrow = length(all_term_ids), ncol = length(curr_logor_list), dimnames = list(all_term_ids, curr_comparisons))
    ora_vi <- matrix(NA, nrow = length(all_term_ids), ncol = length(curr_se_list), dimnames = list(all_term_ids, curr_comparisons))

    for (comp_name in curr_comparisons) {
      ora_yi[names(curr_logor_list[[comp_name]]), comp_name] <- curr_logor_list[[comp_name]]
      ora_vi[names(curr_se_list[[comp_name]]), comp_name] <- curr_se_list[[comp_name]]^2
    }

    curr_op_logor_list <- ora_op_logor_list[[d]]
    curr_op_se_list <- ora_op_se_list[[d]]
    has_operon_data <- operon_enabled && length(curr_op_logor_list) > 0

    ora_op_yi <- matrix(NA, nrow = length(all_term_ids), ncol = length(curr_comparisons), dimnames = list(all_term_ids, curr_comparisons))
    ora_op_vi <- matrix(NA, nrow = length(all_term_ids), ncol = length(curr_comparisons), dimnames = list(all_term_ids, curr_comparisons))

    if (has_operon_data) {
      for (comp_name in names(curr_op_logor_list)) {
        if (comp_name %in% curr_comparisons) {
          ora_op_yi[names(curr_op_logor_list[[comp_name]]), comp_name] <- curr_op_logor_list[[comp_name]]
          ora_op_vi[names(curr_op_se_list[[comp_name]]), comp_name] <- curr_op_se_list[[comp_name]]^2
        }
      }
    }

    cat(sprintf("  Meta-analyzing %d %s across %d comparisons (direction: %s, dual gene+operon)...\n",
                length(all_term_ids), term_noun, length(curr_comparisons), d))

    ora_meta_res <- list()
    for (i in seq_along(all_term_ids)) {
      term_id <- all_term_ids[i]
      yi_term <- ora_yi[i, ]
      vi_term <- ora_vi[i, ]

      valid_term <- !is.na(yi_term) & !is.na(vi_term)
      valid_gse_term <- curr_gse_map[valid_term]
      valid_comp_term <- curr_comparisons[valid_term]
      n_valid_comps <- sum(valid_term)

      if (length(unique(valid_gse_term)) >= max(2, min_ds)) {
        term_sub <- overlaps_by_term[[term_id]]
        if (!is.null(term_sub) && nrow(term_sub) > 0) {
          term_sub <- term_sub[term_sub$comparison %in% valid_comp_term, , drop = FALSE]
        }

        if (!is.null(term_sub) && nrow(term_sub) > 0) {
          gene_comp_counts <- tapply(term_sub$comparison, term_sub$gene, function(x) length(unique(x)))
          sorted_genes <- names(gene_comp_counts)[order(-gene_comp_counts, names(gene_comp_counts))]
          sorted_counts <- as.integer(gene_comp_counts[sorted_genes])

          term_gene_count <- length(sorted_genes)
          term_genes <- paste(sorted_genes, collapse = "/")
          term_genes_symbols <- if (has_symbol_map) paste(.resolve_gene_symbols(sorted_genes, gene_to_symbol), collapse = "/") else NULL
          freq_labels <- .format_gene_freq_labels(sorted_genes, gene_to_symbol)
          term_genes_freq <- paste(sprintf("%s (%d/%d)", freq_labels, sorted_counts, n_valid_comps), collapse = "; ")
        } else {
          term_gene_count <- 0L
          term_genes <- ""
          term_genes_symbols <- if (has_symbol_map) "" else NULL
          term_genes_freq <- ""
        }

        if (!is.null(comparison_control_map)) {
          valid_control_term <- comparison_control_map[valid_comp_term]
          valid_cov_cluster <- paste0(valid_gse_term, "_", ifelse(is.na(valid_control_term), "UnknownControl", valid_control_term))
        } else {
          valid_cov_cluster <- valid_gse_term
        }

        V_mat_logor <- tryCatch(
          metafor::vcalc(vi = vi_term[valid_term], cluster = valid_cov_cluster, obs = 1:sum(valid_term), rho = 0.5),
          error = function(e) vi_term[valid_term]
        )

        if (!debug_printed) {
          cat(sprintf("\n--- DEBUG: Covariance Matrices for %s %s (direction: %s) ---\n", id_col_name, term_id, d))
          cat("Comparisons:\n")
          print(valid_comp_term)
          cat("\nClusters (valid_cov_cluster):\n")
          print(valid_cov_cluster)
          cat("\nCorrelation Matrix (derived from V_mat_logor):\n")
          print(cov2cor(V_mat_logor))
          cat("-----------------------------------------------\n\n")
          debug_printed <- TRUE
        }

        res_re <- tryCatch(
          metafor::rma.mv(yi = yi_term[valid_term], V = V_mat_logor, random = ~ 1 | valid_gse_term / valid_comp_term, method = "REML"),
          error = function(e) {
            tryCatch(
              metafor::rma.mv(yi = yi_term[valid_term], V = V_mat_logor, random = ~ 1 | valid_gse_term / valid_comp_term, method = "REML",
                              control = list(optimizer = "optim", optmethod = "Nelder-Mead")),
              error = function(e2) {
                cat(sprintf("  - [logOR %s] 3-level RE model failed for %s (all optimizers): %s\n", d, term_id, e2$message))
                NULL
              }
            )
          }
        )

        if (!is.null(res_re)) {
          res_re_rob <- tryCatch(
            metafor::robust(res_re, cluster = valid_gse_term, clubSandwich = TRUE),
            error = function(e) {
              cat(sprintf("  - [logOR %s] clubSandwich RVE failed for %s (using model SE): %s\n", d, term_id, e$message))
              NULL
            }
          )
          if (!is.null(res_re_rob)) res_re <- res_re_rob
        }

        res_fe <- tryCatch(
          metafor::rma.mv(yi = yi_term[valid_term], V = V_mat_logor),
          error = function(e) {
            cat(sprintf("  - [logOR FE %s] model failed for %s: %s\n", d, term_id, e$message))
            NULL
          }
        )

        # --- Operon-Level Meta-Analysis (Option B) ---
        op_re_est <- NA_real_
        op_re_se  <- NA_real_
        op_re_p   <- NA_real_
        op_fe_est <- NA_real_
        op_fe_se  <- NA_real_
        op_fe_p   <- NA_real_

        if (has_operon_data && term_id %in% rownames(ora_op_yi)) {
          yi_op <- ora_op_yi[term_id, ]
          vi_op <- ora_op_vi[term_id, ]
          valid_op <- !is.na(yi_op) & !is.na(vi_op)

          if (sum(valid_op) >= 2) {
            valid_gse_op <- curr_gse_map[valid_op]
            valid_comp_op <- curr_comparisons[valid_op]

            if (!is.null(comparison_control_map)) {
              valid_control_op <- comparison_control_map[valid_comp_op]
              valid_cov_cluster_op <- paste0(valid_gse_op, "_", ifelse(is.na(valid_control_op), "UnknownControl", valid_control_op))
            } else {
              valid_cov_cluster_op <- valid_gse_op
            }

            V_mat_op <- tryCatch(
              metafor::vcalc(vi = vi_op[valid_op], cluster = valid_cov_cluster_op, obs = 1:sum(valid_op), rho = 0.5),
              error = function(e) vi_op[valid_op]
            )

            res_re_op <- tryCatch(
              metafor::rma.mv(yi = yi_op[valid_op], V = V_mat_op, random = ~ 1 | valid_gse_op / valid_comp_op, method = "REML"),
              error = function(e) {
                tryCatch(
                  metafor::rma.mv(yi = yi_op[valid_op], V = V_mat_op, random = ~ 1 | valid_gse_op / valid_comp_op, method = "REML",
                                  control = list(optimizer = "optim", optmethod = "Nelder-Mead")),
                  error = function(e2) NULL
                )
              }
            )
            if (!is.null(res_re_op)) {
              res_re_op_rob <- tryCatch(
                metafor::robust(res_re_op, cluster = valid_gse_op, clubSandwich = TRUE),
                error = function(e) NULL
              )
              if (!is.null(res_re_op_rob)) res_re_op <- res_re_op_rob
            }

            res_fe_op <- tryCatch(
              metafor::rma.mv(yi = yi_op[valid_op], V = V_mat_op),
              error = function(e) NULL
            )

            op_re_est <- if (!is.null(res_re_op)) as.numeric(res_re_op$beta)[1] else NA_real_
            op_re_se  <- if (!is.null(res_re_op)) as.numeric(res_re_op$se)[1] else NA_real_
            op_re_p   <- if (!is.null(res_re_op)) .one_sided_p(op_re_est, as.numeric(res_re_op$pval)[1]) else NA_real_

            op_fe_est <- if (!is.null(res_fe_op)) as.numeric(res_fe_op$beta)[1] else NA_real_
            op_fe_se  <- if (!is.null(res_fe_op)) as.numeric(res_fe_op$se)[1] else NA_real_
            op_fe_p   <- if (!is.null(res_fe_op)) .one_sided_p(op_fe_est, as.numeric(res_fe_op$pval)[1]) else NA_real_
          }
        }

        term_op_sub <- op_overlaps_by_term[[term_id]]
        if (!is.null(term_op_sub) && nrow(term_op_sub) > 0) {
          term_op_sub <- term_op_sub[term_op_sub$comparison %in% valid_comp_term, , drop = FALSE]
        }

        if (!is.null(term_op_sub) && nrow(term_op_sub) > 0) {
          op_comp_counts <- tapply(term_op_sub$comparison, term_op_sub$operon_label, function(x) length(unique(x)))
          sorted_ops <- names(op_comp_counts)[order(-op_comp_counts, names(op_comp_counts))]
          sorted_op_counts <- as.integer(op_comp_counts[sorted_ops])

          term_op_count <- length(sorted_ops)
          term_ops <- paste(sorted_ops, collapse = "; ")
          term_ops_freq <- paste(sprintf("%s (%d/%d)", sorted_ops, sorted_op_counts, n_valid_comps), collapse = "; ")
        } else {
          term_op_count <- 0L
          term_ops <- ""
          term_ops_freq <- ""
        }

        if (!is.null(res_re) || !is.null(res_fe)) {
          re_est <- if (!is.null(res_re)) as.numeric(res_re$beta)[1] else NA_real_
          re_se  <- if (!is.null(res_re)) as.numeric(res_re$se)[1] else NA_real_
          re_p   <- if (!is.null(res_re)) .one_sided_p(re_est, as.numeric(res_re$pval)[1]) else NA_real_

          fe_est <- if (!is.null(res_fe)) as.numeric(res_fe$beta)[1] else NA_real_
          fe_se  <- if (!is.null(res_fe)) as.numeric(res_fe$se)[1] else NA_real_
          fe_p   <- if (!is.null(res_fe)) .one_sided_p(fe_est, as.numeric(res_fe$pval)[1]) else NA_real_

          row_data <- data.frame(
            ID = term_id,
            Ontology = term2ontology$ontology[match(term_id, term2ontology$term)],
            Description = term2name$name[match(term_id, term2name$term)],
            Direction = d,
            logOR_RE_estimate = re_est,
            logOR_RE_SE = re_se,
            logOR_RE_pval = re_p,
            logOR_FE_estimate = fe_est,
            logOR_FE_SE = fe_se,
            logOR_FE_pval = fe_p,
            stringsAsFactors = FALSE
          )

          if (has_operon_data) {
            row_data$Operon_logOR_RE_estimate <- op_re_est
            row_data$Operon_logOR_RE_SE <- op_re_se
            row_data$Operon_logOR_RE_pval <- op_re_p
            row_data$Operon_logOR_FE_estimate <- op_fe_est
            row_data$Operon_logOR_FE_SE <- op_fe_se
            row_data$Operon_logOR_FE_pval <- op_fe_p
          }

          row_data$N_comparisons <- n_valid_comps
          row_data$N_datasets <- length(unique(valid_gse_term))
          row_data$Gene_Count <- term_gene_count
          row_data$Genes <- term_genes
          if (has_symbol_map) row_data$Genes_Symbols <- term_genes_symbols
          row_data$Genes_Frequency <- term_genes_freq

          if (has_operon_data) {
            row_data$Operon_Count <- term_op_count
            row_data$Operons <- term_ops
            row_data$Operons_Frequency <- term_ops_freq
          }

          ora_meta_res[[term_id]] <- row_data
        }
      }
    }

    # --- Save and Export Results ---
    if (length(ora_meta_res) > 0) {
      ora_meta_df <- do.call(rbind, ora_meta_res)

      ora_meta_df$logOR_RE_padj <- p.adjust(ora_meta_df$logOR_RE_pval, method = "BH")
      ora_meta_df$logOR_FE_padj <- p.adjust(ora_meta_df$logOR_FE_pval, method = "BH")

      if (has_operon_data) {
        ora_meta_df$Operon_logOR_RE_padj <- p.adjust(ora_meta_df$Operon_logOR_RE_pval, method = "BH")
        ora_meta_df$Operon_logOR_FE_padj <- p.adjust(ora_meta_df$Operon_logOR_FE_pval, method = "BH")

        gene_sig <- !is.na(ora_meta_df$logOR_RE_padj) & ora_meta_df$logOR_RE_padj < report_padj & !is.na(ora_meta_df$logOR_RE_estimate) & ora_meta_df$logOR_RE_estimate > 0
        op_sig   <- !is.na(ora_meta_df$Operon_logOR_RE_padj) & ora_meta_df$Operon_logOR_RE_padj < report_padj & !is.na(ora_meta_df$Operon_logOR_RE_estimate) & ora_meta_df$Operon_logOR_RE_estimate > 0

        ora_meta_df$Operon_Inflation_Flag <- ifelse(
          gene_sig & op_sig,
          "Robust (Gene+Operon Confirmed)",
          ifelse(
            gene_sig & !op_sig,
            ifelse(ora_meta_df$Operon_Count <= 2, "Operon-Driven Inflation (Single/Few Operon)", "Operon-Attenuated (Broadly Dispersed)"),
            ifelse(!gene_sig & op_sig, "Operon-Emergent", "Not Significant")
          )
        )
      }

      out_cols <- c(
        "ID", "Ontology", "Description", "Direction",
        "logOR_RE_estimate", "logOR_RE_SE", "logOR_RE_pval", "logOR_RE_padj",
        "logOR_FE_estimate", "logOR_FE_SE", "logOR_FE_pval", "logOR_FE_padj"
      )
      if (has_operon_data) {
        out_cols <- c(out_cols,
          "Operon_logOR_RE_estimate", "Operon_logOR_RE_SE", "Operon_logOR_RE_pval", "Operon_logOR_RE_padj",
          "Operon_logOR_FE_estimate", "Operon_logOR_FE_SE", "Operon_logOR_FE_pval", "Operon_logOR_FE_padj",
          "Operon_Inflation_Flag"
        )
      }
      out_cols <- c(out_cols,
        "N_comparisons", "N_datasets",
        "Gene_Count", "Genes",
        if (has_symbol_map) "Genes_Symbols",
        "Genes_Frequency"
      )
      if (has_operon_data) {
        out_cols <- c(out_cols, "Operon_Count", "Operons", "Operons_Frequency")
      }
      ora_meta_df <- ora_meta_df[, out_cols, drop = FALSE]
      ora_meta_df <- ora_meta_df[order(ora_meta_df$logOR_RE_pval, na.last = TRUE), ]

      out_file_d <- file.path(out_dir, sprintf("%s_Results_%s.csv", file_prefix, d))
      write.csv(ora_meta_df, out_file_d, row.names = FALSE)
      cat(sprintf("  - Saved %s Results (%s) to %s\n", file_prefix, d, out_file_d))

      if (isTRUE(write_legacy_all) && d == "all") {
        out_file_legacy <- file.path(out_dir, sprintf("%s_Results.csv", file_prefix))
        write.csv(ora_meta_df, out_file_legacy, row.names = FALSE)
        cat(sprintf("  - Saved %s Results (legacy default) to %s\n", file_prefix, out_file_legacy))
      }

      results_by_dir[[d]] <- ora_meta_df

      # --- Detailed Gene Contributions and Gene Pleiotropy Summaries ---
      if (!is.null(all_overlaps_d) && nrow(all_overlaps_d) > 0) {
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

            row_c <- data.frame(
              ID_PLACEHOLDER = t_id,
              Ontology = term2ontology$ontology[match(t_id, term2ontology$term)],
              Description = term2name$name[match(t_id, term2name$term)],
              Direction = d,
              Gene = g_id,
              stringsAsFactors = FALSE
            )
            colnames(row_c)[1] <- id_col_name
            if (has_symbol_map) row_c$Gene_Symbol <- .resolve_gene_symbols(g_id, gene_to_symbol)
            row_c$N_Comparisons <- n_comps
            row_c$Total_Comparisons <- tot_comps
            row_c$Percent_Comparisons <- pct_comps
            row_c$N_Datasets <- length(gses)
            row_c$Comparisons <- paste(sort(comps), collapse = "; ")
            row_c$Datasets <- paste(sort(gses), collapse = "; ")
            row_c
          })
          contrib_df <- do.call(rbind, contrib_rows)
          contrib_df <- contrib_df[order(contrib_df[[id_col_name]], -contrib_df$N_Comparisons, contrib_df$Gene), ]

          contrib_file_d <- file.path(out_dir, sprintf("%s_Gene_Contributions_%s.csv", file_prefix, d))
          write.csv(contrib_df, contrib_file_d, row.names = FALSE)
          cat(sprintf("  - Saved Detailed Gene Contributions (%s) to %s\n", d, contrib_file_d))

          if (isTRUE(write_legacy_all) && d == "all") {
            contrib_file_legacy <- file.path(out_dir, sprintf("%s_Gene_Contributions.csv", file_prefix))
            write.csv(contrib_df, contrib_file_legacy, row.names = FALSE)
            cat(sprintf("  - Saved Detailed Gene Contributions (legacy default) to %s\n", contrib_file_legacy))
          }
          contrib_by_dir[[d]] <- contrib_df

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
              row_s <- data.frame(Gene = g_id, stringsAsFactors = FALSE)
              if (has_symbol_map) row_s$Gene_Symbol <- .resolve_gene_symbols(g_id, gene_to_symbol)
              row_s$Direction <- d
              row_s$N_Significant_Terms <- length(terms)
              row_s$N_Comparisons <- length(comps)
              row_s$N_Datasets <- length(gses)
              row_s$Significant_Terms <- paste(term_desc, collapse = "; ")
              row_s$Comparisons <- paste(sort(comps), collapse = "; ")
              row_s$Datasets <- paste(sort(gses), collapse = "; ")
              row_s
            })
            gene_summary_df <- do.call(rbind, gene_summary_rows)
            gene_summary_df <- gene_summary_df[order(-gene_summary_df$N_Significant_Terms, -gene_summary_df$N_Comparisons, gene_summary_df$Gene), ]

            gene_sum_file_d <- file.path(out_dir, sprintf("%s_Gene_Summary_%s.csv", file_prefix, d))
            write.csv(gene_summary_df, gene_sum_file_d, row.names = FALSE)
            cat(sprintf("  - Saved Gene-Level Summary across significant terms (%s) to %s\n", d, gene_sum_file_d))

            if (isTRUE(write_legacy_all) && d == "all") {
              gene_sum_file_legacy <- file.path(out_dir, sprintf("%s_Gene_Summary.csv", file_prefix))
              write.csv(gene_summary_df, gene_sum_file_legacy, row.names = FALSE)
              cat(sprintf("  - Saved Gene-Level Summary (legacy default) to %s\n", gene_sum_file_legacy))
            }
            summary_by_dir[[d]] <- gene_summary_df
          }
        }
      }

      # --- Detailed Operon Contributions and Operon Summaries ---
      if (has_operon_data && !is.null(all_op_overlaps_d) && nrow(all_op_overlaps_d) > 0) {
        op_overlaps_tested <- all_op_overlaps_d[all_op_overlaps_d$term %in% ora_meta_df$ID, , drop = FALSE]

        if (nrow(op_overlaps_tested) > 0) {
          term_valid_n <- setNames(ora_meta_df$N_comparisons, ora_meta_df$ID)
          to_key <- paste(op_overlaps_tested$term, op_overlaps_tested$operon_id, sep = "@@@")
          to_splits <- split(op_overlaps_tested, to_key)

          op_contrib_rows <- lapply(to_splits, function(df_to) {
            t_id <- df_to$term[1]
            op_id <- df_to$operon_id[1]
            op_lbl <- df_to$operon_label[1]
            comps <- unique(df_to$comparison)
            gses  <- unique(df_to$gse)
            n_comps <- length(comps)
            tot_comps <- if (!is.na(term_valid_n[t_id])) term_valid_n[t_id] else length(curr_comparisons)
            pct_comps <- round(100 * n_comps / tot_comps, 1)

            row_oc <- data.frame(
              ID_PLACEHOLDER = t_id,
              Ontology = term2ontology$ontology[match(t_id, term2ontology$term)],
              Description = term2name$name[match(t_id, term2name$term)],
              Direction = d,
              Operon_ID = op_id,
              Operon_Label = op_lbl,
              N_Comparisons = n_comps,
              Total_Comparisons = tot_comps,
              Percent_Comparisons = pct_comps,
              N_Datasets = length(gses),
              Comparisons = paste(sort(comps), collapse = "; "),
              Datasets = paste(sort(gses), collapse = "; "),
              stringsAsFactors = FALSE
            )
            colnames(row_oc)[1] <- id_col_name
            row_oc
          })
          op_contrib_df <- do.call(rbind, op_contrib_rows)
          op_contrib_df <- op_contrib_df[order(op_contrib_df[[id_col_name]], -op_contrib_df$N_Comparisons, op_contrib_df$Operon_ID), ]

          op_contrib_file_d <- file.path(out_dir, sprintf("%s_Operon_Contributions_%s.csv", file_prefix, d))
          write.csv(op_contrib_df, op_contrib_file_d, row.names = FALSE)
          cat(sprintf("  - Saved Detailed Operon Contributions (%s) to %s\n", d, op_contrib_file_d))

          if (isTRUE(write_legacy_all) && d == "all") {
            op_contrib_file_legacy <- file.path(out_dir, sprintf("%s_Operon_Contributions.csv", file_prefix))
            write.csv(op_contrib_df, op_contrib_file_legacy, row.names = FALSE)
            cat(sprintf("  - Saved Detailed Operon Contributions (legacy default) to %s\n", op_contrib_file_legacy))
          }

          sig_op_term_ids <- ora_meta_df$ID[!is.na(ora_meta_df$Operon_logOR_RE_padj) & ora_meta_df$Operon_logOR_RE_padj < report_padj]
          sig_op_overlaps <- op_overlaps_tested[op_overlaps_tested$term %in% sig_op_term_ids, , drop = FALSE]

          if (nrow(sig_op_overlaps) > 0) {
            op_splits <- split(sig_op_overlaps, sig_op_overlaps$operon_id)
            op_summary_rows <- lapply(op_splits, function(df_o) {
              op_id <- df_o$operon_id[1]
              op_lbl <- df_o$operon_label[1]
              terms <- unique(df_o$term)
              term_names <- term2name$name[match(terms, term2name$term)]
              term_desc <- sprintf("%s (%s)", terms, term_names)
              comps <- unique(df_o$comparison)
              gses <- unique(df_o$gse)
              data.frame(
                Operon_ID = op_id,
                Operon_Label = op_lbl,
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
            op_summary_df <- do.call(rbind, op_summary_rows)
            op_summary_df <- op_summary_df[order(-op_summary_df$N_Significant_Terms, -op_summary_df$N_Comparisons, op_summary_df$Operon_ID), ]

            op_sum_file_d <- file.path(out_dir, sprintf("%s_Operon_Summary_%s.csv", file_prefix, d))
            write.csv(op_summary_df, op_sum_file_d, row.names = FALSE)
            cat(sprintf("  - Saved Operon-Level Summary across significant terms (%s) to %s\n", d, op_sum_file_d))

            if (isTRUE(write_legacy_all) && d == "all") {
              op_sum_file_legacy <- file.path(out_dir, sprintf("%s_Operon_Summary.csv", file_prefix))
              write.csv(op_summary_df, op_sum_file_legacy, row.names = FALSE)
              cat(sprintf("  - Saved Operon-Level Summary (legacy default) to %s\n", op_sum_file_legacy))
            }
          }
        }
      }

      # --- Reporting-only semantic clustering (GO only) ---
      if (isTRUE(run_clustering)) {
        clustered <- tryCatch(
          cluster_significant_go_terms(ora_meta_df, sig_padj = report_padj, sim_cutoff = sem_sim_cutoff),
          error = function(e) {
            cat(sprintf("  - [%s] GOSemSim clustering failed: %s\n", d, e$message))
            NULL
          }
        )
        if (!is.null(clustered)) {
          clustered <- clustered[order(clustered$Ontology, clustered$Cluster, !clustered$Is_Representative, clustered$logOR_RE_padj), ]
          clusters_file <- file.path(out_dir, sprintf("%s_Significant_Term_Clusters_%s.csv", file_prefix, d))
          write.csv(clustered, clusters_file, row.names = FALSE)
          reps <- clustered[clustered$Is_Representative, ]
          reps <- reps[order(reps$logOR_RE_padj), ]
          reps_file <- file.path(out_dir, sprintf("%s_Representative_Terms_%s.csv", file_prefix, d))
          write.csv(reps, reps_file, row.names = FALSE)
          cat(sprintf("  - %d significant terms (padj < %g) -> %d clusters (Wang sim >= %g). Saved clusters to %s and representatives to %s\n",
                      nrow(clustered), report_padj, nrow(reps), sem_sim_cutoff, clusters_file, reps_file))
        } else {
          cat(sprintf("  - [%s] No significant terms (padj < %g) to cluster.\n", d, report_padj))
        }
      }
    } else {
      cat(sprintf("  - [%s] Not enough valid datasets/comparisons to run meta-analysis on any %s.\n", d, term_noun))
    }
  }

  if (length(results_by_dir) > 1) {
    combined_df <- do.call(rbind, results_by_dir)
    combined_file <- file.path(out_dir, sprintf("%s_Results_combined.csv", file_prefix))
    write.csv(combined_df, combined_file, row.names = FALSE)
    cat(sprintf("  - Saved Combined %s Results across all directions to %s\n", file_prefix, combined_file))

    if (length(contrib_by_dir) > 0) {
      combined_contrib <- do.call(rbind, contrib_by_dir)
      combined_contrib_file <- file.path(out_dir, sprintf("%s_Gene_Contributions_combined.csv", file_prefix))
      write.csv(combined_contrib, combined_contrib_file, row.names = FALSE)
      cat(sprintf("  - Saved Combined Detailed Gene Contributions across all directions to %s\n", combined_contrib_file))
    }

    if (length(summary_by_dir) > 0) {
      combined_summary <- do.call(rbind, summary_by_dir)
      combined_sum_file <- file.path(out_dir, sprintf("%s_Gene_Summary_combined.csv", file_prefix))
      write.csv(combined_summary, combined_sum_file, row.names = FALSE)
      cat(sprintf("  - Saved Combined Gene-Level Summary across all directions to %s\n", combined_sum_file))
    }
  }

  if (length(results_by_dir) == 0) {
    cat(sprintf("  - No %s meta-analysis results could be computed for any direction.\n", file_prefix))
  }

  invisible(results_by_dir)
}

# Shared Fixed-Effect (FE) significant-gene ORA engine + side-by-side comparison with Meta-ORA.
.run_geneset_fe_ora_core <- function(prep, res_df, out_dir, ora_padj = 0.05, ora_lfc = 0.5, directions = c("all", "up", "down"), sem_sim_cutoff = 0.7, report_padj = 0.05, operon_mapping_file = "operon_gene_mapping.csv", operon_de_threshold = 0.5, meta_results_by_dir = NULL, file_prefix = "FE_ORA", cmp_filename = "FE_vs_Meta_ORA_Comparison.csv", term_noun = "terms", run_clustering = TRUE) {
  t2g_u <- prep$term2gene
  term2name <- prep$term2name
  term2ontology <- prep$term2ontology
  valid_terms <- prep$valid_terms
  cluster_id_to_symbol <- prep$cluster_id_to_gene
  gene_to_symbol <- prep$gene_to_symbol
  has_symbol_map <- !is.null(gene_to_symbol)

  if (length(valid_terms) == 0 || length(prep$global_genes) == 0) return(invisible(NULL))

  universe_clusters <- as.character(res_df$cluster_id)
  N <- length(prep$global_genes)
  K <- as.numeric(table(factor(t2g_u$term, levels = valid_terms)))

  op_info <- .setup_operon_units(operon_mapping_file, universe_clusters, cluster_id_to_symbol, t2g_u)
  operon_enabled <- op_info$enabled
  if (operon_enabled) {
    op_to_cids <- op_info$op_to_cids
    op_labels <- op_info$op_labels
    term_to_ops <- op_info$term_to_ops
    op_meas <- sapply(op_to_cids, function(cids) sum(cids %in% universe_clusters))
    valid_comp_ops <- names(op_to_cids)[op_meas > 0]
    bg_op_num <- sapply(valid_terms, function(t) { o <- term_to_ops[[t]]; if (is.null(o)) 0L else sum(valid_comp_ops %in% o) })
  }

  hg_table <- function(a, de_n, K_vec, N_tot) {
    b <- de_n - a; cc <- K_vec - a; dd <- N_tot - de_n - cc
    p <- phyper(a - 1, K_vec, N_tot - K_vec, de_n, lower.tail = FALSE)
    lor <- log(((a + 0.5) * (dd + 0.5)) / ((b + 0.5) * (cc + 0.5)))
    se <- sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (cc + 0.5) + 1 / (dd + 0.5))
    list(p = p, lor = lor, se = se, fold = (a / de_n) / (K_vec / N_tot))
  }

  results_by_dir <- list()
  for (d in directions) {
    sig <- res_df$Meta_FE_padj < ora_padj & switch(d,
      all = abs(res_df$Meta_FE_LFC) > ora_lfc,
      up = res_df$Meta_FE_LFC > ora_lfc,
      down = res_df$Meta_FE_LFC < -ora_lfc)
    sig[is.na(sig)] <- FALSE
    de_clusters <- as.character(res_df$cluster_id[sig])
    de_symbols <- unique(cluster_id_to_symbol[de_clusters])
    de_symbols <- de_symbols[!is.na(de_symbols) & de_symbols != ""]
    de_n <- length(de_symbols)
    cat(sprintf("  [%s] FE-significant clusters: %d | mapped DE genes: %d\n", d, length(de_clusters), de_n))
    if (de_n == 0) next

    ov <- t2g_u[t2g_u$gene %in% de_symbols, , drop = FALSE]
    a <- as.numeric(table(factor(ov$term, levels = valid_terms)))
    keep <- K > 0
    h <- hg_table(a, de_n, K, N)
    df <- data.frame(
      ID = valid_terms, Ontology = term2ontology$ontology[match(valid_terms, term2ontology$term)],
      Description = term2name$name[match(valid_terms, term2name$term)], Direction = d,
      Gene_Count = a, Term_Size = K, DE_Size = de_n, Universe_Size = N,
      Fold_Enrichment = h$fold, logOR_estimate = h$lor, logOR_SE = h$se, pval = h$p,
      stringsAsFactors = FALSE
    )[keep, , drop = FALSE]
    df$padj <- p.adjust(df$pval, method = "BH")
    ov_genes <- split(ov$gene, ov$term)
    df$Genes <- sapply(df$ID, function(t) paste(sort(unique(ov_genes[[t]])), collapse = "/"))
    if (has_symbol_map) {
      df$Genes_Symbols <- sapply(df$ID, function(t) {
        gs <- sort(unique(ov_genes[[t]]))
        paste(.resolve_gene_symbols(gs, gene_to_symbol), collapse = "/")
      })
    }

    if (operon_enabled) {
      de_ops_all <- valid_comp_ops[sapply(valid_comp_ops, function(o) sum(op_to_cids[[o]] %in% de_clusters) / op_meas[[o]]) >= operon_de_threshold]
      de_op_n <- length(de_ops_all); tot_ops <- length(valid_comp_ops)
      if (de_op_n > 0) {
        a_op <- sapply(valid_terms, function(t) { o <- term_to_ops[[t]]; if (is.null(o)) 0L else sum(de_ops_all %in% o) })
        ho <- hg_table(a_op, de_op_n, bg_op_num, tot_ops)
        df$Operon_Count <- a_op[keep]; df$Operon_Term_Size <- bg_op_num[keep]
        df$Operon_logOR_estimate <- ho$lor[keep]; df$Operon_logOR_SE <- ho$se[keep]
        df$Operon_pval <- ho$p[keep]; df$Operon_padj <- p.adjust(df$Operon_pval, method = "BH")
        df$Operons <- sapply(df$ID, function(t) { o <- intersect(de_ops_all, term_to_ops[[t]]); paste(unname(op_labels[o]), collapse = "; ") })
        gs <- !is.na(df$padj) & df$padj < report_padj & df$logOR_estimate > 0
        os <- !is.na(df$Operon_padj) & df$Operon_padj < report_padj & df$Operon_logOR_estimate > 0
        df$Operon_Inflation_Flag <- ifelse(gs & os, "Robust (Gene+Operon Confirmed)",
          ifelse(gs & !os, ifelse(df$Operon_Count <= 2, "Operon-Driven Inflation (Single/Few Operon)", "Operon-Attenuated (Broadly Dispersed)"),
            ifelse(!gs & os, "Operon-Emergent", "Not Significant")))
      }
    }
    df <- df[order(df$pval), ]
    f <- file.path(out_dir, sprintf("%s_Results_%s.csv", file_prefix, d))
    write.csv(df, f, row.names = FALSE)
    cat(sprintf("  - Saved %s results (%s): %d %s tested, %d significant (padj < %g) -> %s\n",
                file_prefix, d, nrow(df), term_noun, sum(df$padj < report_padj & df$logOR_estimate > 0, na.rm = TRUE), report_padj, f))
    results_by_dir[[d]] <- df

    if (isTRUE(run_clustering)) {
      tmp <- df; tmp$logOR_RE_padj <- ifelse(tmp$logOR_estimate > 0, tmp$padj, NA_real_); tmp$logOR_RE_estimate <- tmp$logOR_estimate
      clustered <- tryCatch(cluster_significant_go_terms(tmp, sig_padj = report_padj, sim_cutoff = sem_sim_cutoff),
        error = function(e) { cat(sprintf("  - [%s] GOSemSim clustering failed: %s\n", d, e$message)); NULL })
      if (!is.null(clustered)) {
        clustered <- clustered[order(clustered$Ontology, clustered$Cluster, !clustered$Is_Representative, clustered$padj), ]
        write.csv(clustered, file.path(out_dir, sprintf("%s_Significant_Term_Clusters_%s.csv", file_prefix, d)), row.names = FALSE)
        write.csv(clustered[clustered$Is_Representative, ][order(clustered$padj[clustered$Is_Representative]), ],
          file.path(out_dir, sprintf("%s_Representative_Terms_%s.csv", file_prefix, d)), row.names = FALSE)
      }
    }
  }

  # --- Side-by-side comparison with meta-ORA ---
  if (!is.null(meta_results_by_dir) && length(results_by_dir) > 0) {
    cmp <- list()
    for (d in names(results_by_dir)) {
      fe <- results_by_dir[[d]]; me <- meta_results_by_dir[[d]]
      if (is.null(me)) next
      all_ids <- union(fe$ID, me$ID)
      m <- data.frame(ID = all_ids, Direction = d, stringsAsFactors = FALSE)
      m$Description <- c(fe$Description, me$Description)[match(all_ids, c(fe$ID, me$ID))]
      m$Ontology <- c(fe$Ontology, me$Ontology)[match(all_ids, c(fe$ID, me$ID))]
      i <- match(all_ids, fe$ID); j <- match(all_ids, me$ID)
      m$FE_ORA_logOR <- fe$logOR_estimate[i]; m$FE_ORA_padj <- fe$padj[i]
      m$MetaORA_RE_logOR <- me$logOR_RE_estimate[j]; m$MetaORA_RE_padj <- me$logOR_RE_padj[j]
      m$MetaORA_FE_logOR <- me$logOR_FE_estimate[j]; m$MetaORA_FE_padj <- me$logOR_FE_padj[j]
      sf <- !is.na(m$FE_ORA_padj) & m$FE_ORA_padj < report_padj & m$FE_ORA_logOR > 0
      sm <- !is.na(m$MetaORA_RE_padj) & m$MetaORA_RE_padj < report_padj & m$MetaORA_RE_logOR > 0
      m$Concordance <- ifelse(sf & sm, "Both", ifelse(sf, "FE_ORA_only", ifelse(sm, "MetaORA_only", "Neither")))
      cmp[[d]] <- m[order(m$Concordance, m$FE_ORA_padj), ]
      cat(sprintf("  [%s] significant %s: Both=%d | FE_ORA only=%d | Meta-ORA only=%d\n", d, term_noun, sum(sf & sm), sum(sf & !sm), sum(sm & !sf)))
    }
    if (length(cmp) > 0) {
      cmp_df <- do.call(rbind, cmp)
      write.csv(cmp_df, file.path(out_dir, cmp_filename), row.names = FALSE)
      cat(sprintf("  - Saved %s\n", cmp_filename))
    }
  }
  invisible(results_by_dir)
}

run_ora_meta_analysis <- function(gaf_file, lfc_list, pval_list, se_list, cts, cluster_annotation, comparison_gse_map, comparison_control_map = NULL, min_ds, out_dir, ora_padj = 0.05, ora_lfc = 1.0, exclude_evidence = c("ND"), min_gs_size = 5, max_gs_size = 500, directions = c("all", "up", "down"), sem_sim_cutoff = 0.7, report_padj = 0.05, operon_mapping_file = "operon_gene_mapping.csv", operon_de_threshold = 0.5) {
  cat("\n=== RUNNING OVER-REPRESENTATION META-ANALYSIS ===\n")
  directions <- unique(match.arg(directions, c("all", "up", "down"), several.ok = TRUE))
  cat(sprintf("Directions to analyze: %s\n", paste(directions, collapse = ", ")))

  if (is.null(gaf_file) || !file.exists(gaf_file)) {
    cat("  - GAF file not found or not provided. Skipping ORA meta-analysis.\n")
    return(invisible(NULL))
  }

  cat("Loading GAF file:", gaf_file, "\n")
  prep <- .prepare_go_annotations(
    gaf_file = gaf_file,
    cluster_annotation = cluster_annotation,
    study_clusters = rownames(cts),
    exclude_evidence = exclude_evidence,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size
  )

  .run_geneset_meta_ora_core(
    prep = prep,
    lfc_list = lfc_list,
    pval_list = pval_list,
    se_list = se_list,
    cts = cts,
    comparison_gse_map = comparison_gse_map,
    comparison_control_map = comparison_control_map,
    min_ds = min_ds,
    out_dir = out_dir,
    ora_padj = ora_padj,
    ora_lfc = ora_lfc,
    directions = directions,
    sem_sim_cutoff = sem_sim_cutoff,
    report_padj = report_padj,
    operon_mapping_file = operon_mapping_file,
    operon_de_threshold = operon_de_threshold,
    file_prefix = "ORA_Meta_Analysis",
    id_col_name = "GO_ID",
    term_noun = "GO terms",
    write_legacy_all = TRUE,
    run_clustering = TRUE
  )
}

# ======================================================================================
# ORA on the Fixed-Effect (FE) meta-analysis significant gene lists
# ======================================================================================
run_fe_ora <- function(res_df, gaf_file, cluster_annotation, out_dir, ora_padj = 0.05, ora_lfc = 0.5, exclude_evidence = c("ND"), min_gs_size = 5, max_gs_size = 500, directions = c("all", "up", "down"), sem_sim_cutoff = 0.7, report_padj = 0.05, operon_mapping_file = "operon_gene_mapping.csv", operon_de_threshold = 0.5, meta_results_by_dir = NULL) {
  cat("\n=== RUNNING ORA ON FE-SIGNIFICANT GENES ===\n")
  if (is.null(gaf_file) || !file.exists(gaf_file)) {
    cat("  - GAF file not found or not provided. Skipping FE ORA.\n")
    return(invisible(NULL))
  }
  directions <- unique(match.arg(directions, c("all", "up", "down"), several.ok = TRUE))

  if (!"cluster_id" %in% colnames(res_df)) res_df$cluster_id <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$Meta_FE_pval) & !is.na(res_df$Meta_FE_LFC), , drop = FALSE]
  universe_clusters <- as.character(res_df$cluster_id)

  prep <- .prepare_go_annotations(
    gaf_file = gaf_file,
    cluster_annotation = cluster_annotation,
    study_clusters = universe_clusters,
    exclude_evidence = exclude_evidence,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size
  )
  cat(sprintf("  FE universe: %d clusters -> %d mapped GAF genes\n", length(universe_clusters), length(prep$global_genes)))

  .run_geneset_fe_ora_core(
    prep = prep,
    res_df = res_df,
    out_dir = out_dir,
    ora_padj = ora_padj,
    ora_lfc = ora_lfc,
    directions = directions,
    sem_sim_cutoff = sem_sim_cutoff,
    report_padj = report_padj,
    operon_mapping_file = operon_mapping_file,
    operon_de_threshold = operon_de_threshold,
    meta_results_by_dir = meta_results_by_dir,
    file_prefix = "FE_ORA",
    cmp_filename = "FE_vs_Meta_ORA_Comparison.csv",
    term_noun = "terms",
    run_clustering = TRUE
  )
}

# ======================================================================================
# KEGG Pathway Meta-ORA and FE-ORA (S. aureus NCTC 8325 = "sao")
# ======================================================================================
run_kegg_meta_analysis <- function(lfc_list, pval_list, se_list, cts, cluster_annotation, comparison_gse_map, comparison_control_map = NULL, min_ds, out_dir, ora_padj = 0.05, ora_lfc = 1.0, organism = "sao", kegg_cache_file = "kegg_sao_cache.rds", refresh_kegg = FALSE, min_gs_size = 5, max_gs_size = 500, directions = c("all", "up", "down"), report_padj = 0.05, operon_mapping_file = "operon_gene_mapping.csv", operon_de_threshold = 0.5, exclude_overview_011 = TRUE) {
  cat("\n=== RUNNING KEGG PATHWAY OVER-REPRESENTATION META-ANALYSIS ===\n")
  directions <- unique(match.arg(directions, c("all", "up", "down"), several.ok = TRUE))
  cat(sprintf("Directions to analyze: %s\n", paste(directions, collapse = ", ")))

  prep <- .prepare_kegg_annotations(
    cluster_annotation = cluster_annotation,
    study_clusters = rownames(cts),
    organism = organism,
    kegg_cache_file = kegg_cache_file,
    refresh_kegg = refresh_kegg,
    operon_mapping_file = operon_mapping_file,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size,
    exclude_overview_011 = exclude_overview_011
  )

  .run_geneset_meta_ora_core(
    prep = prep,
    lfc_list = lfc_list,
    pval_list = pval_list,
    se_list = se_list,
    cts = cts,
    comparison_gse_map = comparison_gse_map,
    comparison_control_map = comparison_control_map,
    min_ds = min_ds,
    out_dir = out_dir,
    ora_padj = ora_padj,
    ora_lfc = ora_lfc,
    directions = directions,
    sem_sim_cutoff = 0.7,
    report_padj = report_padj,
    operon_mapping_file = operon_mapping_file,
    operon_de_threshold = operon_de_threshold,
    file_prefix = "KEGG_ORA_Meta_Analysis",
    id_col_name = "KEGG_ID",
    term_noun = "KEGG pathways",
    write_legacy_all = FALSE,
    run_clustering = FALSE
  )
}

run_fe_kegg_ora <- function(res_df, cluster_annotation, out_dir, ora_padj = 0.05, ora_lfc = 0.5, organism = "sao", kegg_cache_file = "kegg_sao_cache.rds", refresh_kegg = FALSE, min_gs_size = 5, max_gs_size = 500, directions = c("all", "up", "down"), report_padj = 0.05, operon_mapping_file = "operon_gene_mapping.csv", operon_de_threshold = 0.5, meta_results_by_dir = NULL, exclude_overview_011 = TRUE) {
  cat("\n=== RUNNING KEGG PATHWAY ORA ON FE-SIGNIFICANT GENES ===\n")
  directions <- unique(match.arg(directions, c("all", "up", "down"), several.ok = TRUE))

  if (!"cluster_id" %in% colnames(res_df)) res_df$cluster_id <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$Meta_FE_pval) & !is.na(res_df$Meta_FE_LFC), , drop = FALSE]
  universe_clusters <- as.character(res_df$cluster_id)

  prep <- .prepare_kegg_annotations(
    cluster_annotation = cluster_annotation,
    study_clusters = universe_clusters,
    organism = organism,
    kegg_cache_file = kegg_cache_file,
    refresh_kegg = refresh_kegg,
    operon_mapping_file = operon_mapping_file,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size,
    exclude_overview_011 = exclude_overview_011
  )
  cat(sprintf("  FE universe: %d clusters -> %d mapped KEGG genes\n", length(universe_clusters), length(prep$global_genes)))

  .run_geneset_fe_ora_core(
    prep = prep,
    res_df = res_df,
    out_dir = out_dir,
    ora_padj = ora_padj,
    ora_lfc = ora_lfc,
    directions = directions,
    sem_sim_cutoff = 0.7,
    report_padj = report_padj,
    operon_mapping_file = operon_mapping_file,
    operon_de_threshold = operon_de_threshold,
    meta_results_by_dir = meta_results_by_dir,
    file_prefix = "FE_KEGG_ORA",
    cmp_filename = "FE_vs_Meta_KEGG_ORA_Comparison.csv",
    term_noun = "KEGG pathways",
    run_clustering = FALSE
  )
}
