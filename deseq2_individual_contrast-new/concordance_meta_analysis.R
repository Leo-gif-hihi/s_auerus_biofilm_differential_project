# ============================================================================
# Consistently CONCORDANT / DISCORDANT gene groups across datasets
# (random-effects meta-analysis of gene-gene association)
# ============================================================================
#
# HOW TO RUN
#   source("concordance_meta_analysis.R")
#   run_concordance_meta_analysis(
#     dataset_dds_list   = dataset_dds_list,
#     out_dir            = out_dir,
#     cluster_annotation = cluster_annotation,
#     mode               = "with_group"  # or "within_group"
#   )
#
# IDEA
#   For every gene pair (A, B) and every dataset k:
#     r_k  = correlation of A and B across the samples of dataset k
#            (VST expression, batch removed)
#     z_k  = atanh(r_k),  var(z_k) = 1 / (n_k - p - 3)        [Fisher z]
#   Then pool z_k across datasets with a RANDOM-EFFECTS model (DerSimonian-Laird
#   tau^2, Hartung-Knapp inference):
#     pooled r > 0, low heterogeneity, same sign in every dataset -> CONCORDANT
#     pooled r < 0, low heterogeneity, same sign in every dataset -> DISCORDANT
#   Heterogeneity (tau^2, I^2) and sign counts are what make it "consistent".
#
#   mode = "with_group"   : keep the biofilm/control shift in the data, so
#                           r_k reflects whether A and B MOVE TOGETHER between
#                           conditions (closest to "same direction of DE").
#   mode = "within_group" : remove the group effect first -> co-regulation
#                           independent of the treatment (needs larger n).
#
#   NOTE: the sign multiplier you apply to LFCs in the main script is irrelevant
#   here. Correlations do not depend on the contrast orientation.
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(metafor)
})

fisher_z <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))

# ---- 1. per-dataset expression (VST, batch removed) ------------------------
prep_expr <- function(dds, mode = c("with_group", "within_group")) {
  mode <- match.arg(mode)
  cd   <- as.data.frame(colData(dds))

  # Robust VST estimation with graceful fallback for small samples or low counts
  vsd <- tryCatch(
    assay(vst(dds, blind = FALSE)),
    error = function(e) {
      tryCatch(
        assay(varianceStabilizingTransformation(dds, blind = FALSE)),
        error = function(e2) {
          # Fallback when sample size or dispersion fit prevents VST
          log2(counts(dds, normalized = TRUE) + 1)
        }
      )
    }
  )

  batch <- droplevels(factor(cd$Batch))
  group <- droplevels(factor(cd$Group))

  if (mode == "with_group") {
    # remove batch, protect the group effect
    x <- if (nlevels(batch) > 1) {
      tryCatch(
        removeBatchEffect(vsd, batch = batch, design = model.matrix(~group)),
        error = function(e) vsd
      )
    } else {
      vsd
    }
    p <- max(0, nlevels(batch) - 1)
  } else {
    # remove batch AND group (residuals)
    des <- if (nlevels(batch) > 1) model.matrix(~ batch + group) else model.matrix(~group)
    x   <- tryCatch(
      residuals(lmFit(vsd, des), vsd),
      error = function(e) vsd
    )
    p   <- ncol(des) - 1
  }

  # df = n - p - 3 is the inverse of Var(Fisher z) for a partial correlation
  list(x = x, df = ncol(x) - p - 3)
}

# ---- 2. all-pairs random-effects meta-analysis (vectorised, 2 passes) ------
coexpr_meta <- function(prep_list, min_k = 3, use_hk = TRUE) {
  prep_list <- prep_list[sapply(prep_list, function(p) !is.null(p$df) && p$df >= 1)]
  if (length(prep_list) < min_k) {
    cat(sprintf("  [Warning] Fewer than min_k (%d) datasets have df >= 1 (found %d). Cannot perform concordance meta-analysis.\n",
                min_k, length(prep_list)))
    return(NULL)
  }

  genes <- sort(unique(unlist(lapply(prep_list, function(p) rownames(p$x)))))
  G     <- length(genes)
  cat(sprintf("  Genes: %d | Datasets used with df >= 1: %d\n", G, length(prep_list)))
  zero  <- function() matrix(0, G, G, dimnames = list(genes, genes))

  zmat_for <- function(p) {
    x   <- p$x
    ok  <- apply(x, 1, sd) > 0
    ok[is.na(ok)] <- FALSE
    r   <- cor(t(x[ok, , drop = FALSE]))
    Z   <- matrix(NA_real_, G, G, dimnames = list(genes, genes))
    idx <- match(rownames(r), genes)
    Z[idx, idx] <- fisher_z(r)
    Z
  }

  # pass 1: fixed-effect sums -> Q, tau^2 (DerSimonian-Laird)
  K <- zero(); W <- zero(); W2 <- zero(); WZ <- zero(); WZ2 <- zero()
  NP <- zero(); NN <- zero()
  for (p in prep_list) {
    Z <- zmat_for(p); ok <- !is.na(Z); Z[!ok] <- 0; w <- p$df
    K   <- K   + ok
    W   <- W   + w * ok
    W2  <- W2  + w^2 * ok
    WZ  <- WZ  + w * Z
    WZ2 <- WZ2 + w * Z^2
    NP  <- NP  + (Z > 0)
    NN  <- NN  + (Z < 0)
  }
  Q    <- WZ2 - WZ^2 / W
  Cc   <- W - W2 / W
  tau2 <- pmax(0, (Q - (K - 1)) / Cc)
  tau2[!is.finite(tau2)] <- 0
  rm(W2, WZ, WZ2, Cc); gc()

  # pass 2: random-effects weights
  Ws <- zero(); WsZ <- zero(); WsZ2 <- zero()
  for (p in prep_list) {
    Z <- zmat_for(p); ok <- !is.na(Z); Z[!ok] <- 0
    ws   <- ok / (1 / p$df + tau2)
    Ws   <- Ws   + ws
    WsZ  <- WsZ  + ws * Z
    WsZ2 <- WsZ2 + ws * Z^2
  }
  zhat <- WsZ / Ws
  if (use_hk) {
    qhk <- (WsZ2 - WsZ^2 / Ws) / pmax(K - 1, 1)
    qhk[!is.finite(qhk)] <- 1
    se  <- sqrt(pmax(1, qhk) / Ws)           # Hartung-Knapp (floored at 1)
    pv  <- zero(); pv[] <- 2 * pt(-abs(zhat / se), df = pmax(K - 1, 1))
  } else {
    se  <- 1 / sqrt(Ws)
    pv  <- zero(); pv[] <- 2 * pnorm(-abs(zhat / se))
  }
  I2 <- pmax(0, (Q - (K - 1)) / Q) * 100
  I2[!is.finite(I2)] <- 0

  list(r = tanh(zhat), z = zhat, se = se, p = pv, I2 = I2, tau2 = tau2,
       K = K, NP = NP, NN = NN)
}

# ---- 3. call consistent concordant / discordant gene pairs -----------------
call_edges <- function(res, min_k = 3, fdr = 0.05, min_abs_r = 0.5, max_I2 = 50) {
  if (is.null(res)) return(data.frame())
  sel <- upper.tri(res$r) & res$K >= min_k & !is.na(res$p)
  idx <- which(sel)
  if (length(idx) == 0) return(data.frame())

  padj <- p.adjust(res$p[idx], method = "BH")        # BH over ALL tested pairs
  same_sign <- res$NP[idx] == 0 | res$NN[idx] == 0   # same direction in every dataset
  pass <- padj < fdr & abs(res$r[idx]) >= min_abs_r & res$I2[idx] <= max_I2 & same_sign
  padj <- padj[pass]; idx <- idx[pass]
  if (length(idx) == 0) return(data.frame())

  ij <- arrayInd(idx, dim(res$r))
  data.frame(
    gene1 = rownames(res$r)[ij[, 1]],
    gene2 = colnames(res$r)[ij[, 2]],
    k_datasets = res$K[idx],
    pooled_r   = res$r[idx],
    se_z       = res$se[idx],
    padj       = padj,
    I2         = res$I2[idx],
    tau2       = res$tau2[idx],
    n_datasets_pos = res$NP[idx],
    n_datasets_neg = res$NN[idx],
    relationship = ifelse(res$r[idx] > 0, "concordant", "discordant"),
    stringsAsFactors = FALSE
  )
}

# ---- 4. groups of genes: modules from consistent concordant structure ------
find_modules <- function(res, min_k = 3, max_I2 = 50, cut_h = 0.6, min_size = 5) {
  if (is.null(res)) return(character(0))
  R <- res$r
  bad <- is.na(R) | res$K < min_k | res$I2 > max_I2 | (res$NP > 0 & res$NN > 0)
  R[bad] <- 0
  diag(R) <- 1
  keep <- rowSums(abs(R) >= 0.5) > 1
  if (sum(keep) < max(2, min_size)) {
    cat(sprintf("  Fewer than %d genes passed module connectivity criteria. Skipping module clustering.\n",
                max(2, min_size)))
    return(character(0))
  }
  R <- R[keep, keep, drop = FALSE]
  hc <- hclust(as.dist(1 - R), method = "average")
  cl <- cutree(hc, h = cut_h)
  tab <- table(cl)
  good <- as.integer(names(tab)[tab >= min_size])
  if (length(good) == 0) {
    cat(sprintf("  No modules met minimum size threshold (min_size = %d).\n", min_size))
    return(character(0))
  }
  cl[!cl %in% good] <- NA
  lab <- setNames(paste0("M", match(cl, good)), names(cl))
  lab[is.na(cl)] <- NA
  lab <- lab[!is.na(lab)]
  lab
}

# ---- 5. module-vs-module random-effects test (few tests, less noise) -------
module_meta <- function(prep_list, modules, max_I2 = 50) {
  mods <- sort(unique(na.omit(modules)))
  if (length(mods) < 2) return(NULL)
  per_ds <- lapply(prep_list, function(p) {
    if (is.null(p) || p$df < 1) return(NULL)
    x  <- t(scale(t(p$x)))
    sc <- sapply(mods, function(m) {
      g <- intersect(names(modules)[which(modules == m)], rownames(x))
      if (length(g) < 3) rep(NA_real_, ncol(x)) else colMeans(x[g, , drop = FALSE], na.rm = TRUE)
    })
    colnames(sc) <- mods
    list(cor = cor(sc, use = "pairwise.complete.obs"), df = p$df)
  })
  out <- list()
  for (pr in combn(length(mods), 2, simplify = FALSE)) {
    i <- pr[1]; j <- pr[2]
    z <- sapply(per_ds, function(d) if (is.null(d)) NA_real_ else fisher_z(d$cor[i, j]))
    v <- sapply(per_ds, function(d) if (is.null(d)) NA_real_ else 1 / d$df)
    ok <- !is.na(z) & !is.na(v)
    if (sum(ok) < 3) next
    fit <- tryCatch(rma(yi = z[ok], vi = v[ok], method = "REML", test = "knha"),
                    error = function(e) NULL)
    if (is.null(fit)) next
    pi <- tryCatch(predict(fit), error = function(e) NULL)
    out[[length(out) + 1]] <- data.frame(
      module1 = mods[i], module2 = mods[j], k_datasets = sum(ok),
      pooled_r = tanh(as.numeric(fit$beta)),
      ci_lb = tanh(fit$ci.lb), ci_ub = tanh(fit$ci.ub),
      pi_lb = if (!is.null(pi)) tanh(pi$pi.lb) else NA_real_,
      pi_ub = if (!is.null(pi)) tanh(pi$pi.ub) else NA_real_,
      pval = fit$pval, I2 = fit$I2, tau2 = fit$tau2,
      n_datasets_pos = sum(z[ok] > 0), n_datasets_neg = sum(z[ok] < 0),
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(NULL)
  df <- do.call(rbind, out)
  df$padj <- p.adjust(df$pval, "BH")
  df$relationship <- ifelse(df$pooled_r > 0, "concordant", "discordant")
  df$consistent <- (df$n_datasets_pos == 0 | df$n_datasets_neg == 0) & df$I2 <= max_I2
  df
}

# ============================================================================
# MAIN PIPELINE ENTRY POINT
# ============================================================================
run_concordance_meta_analysis <- function(
  dataset_dds_list,
  out_dir            = ".",
  cluster_annotation = NULL,
  mode               = "with_group",
  min_k              = 3,
  fdr_cut            = 0.05,
  min_abs_r          = 0.5,
  max_I2             = 50,
  module_cut_h       = 0.6,
  module_min_sz      = 5
) {
  cat("\n" , paste(rep("=", 70), collapse = ""), "\n")
  cat("GENE-PAIR CONCORDANCE META-ANALYSIS (mode:", mode, ")\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")

  if (is.null(dataset_dds_list) || length(dataset_dds_list) == 0) {
    cat("  [Notice] No valid DESeqDataSet objects found in dataset_dds_list. Skipping concordance meta-analysis.\n")
    return(NULL)
  }

  cat("  Pre-processing datasets (VST + batch removal)...\n")
  prep_list <- lapply(dataset_dds_list, prep_expr, mode = mode)
  df_summary <- sapply(prep_list, function(p) if (is.null(p)) NA else p$df)
  cat("  Effective degrees of freedom (df = n - p - 3) per dataset:\n")
  print(df_summary)

  res <- coexpr_meta(prep_list, min_k = min_k)
  if (is.null(res)) {
    cat("  Co-expression meta-analysis could not be completed (insufficient datasets with df >= 1).\n")
    return(NULL)
  }

  # --- 1. Call consistent pairs ---
  edges <- call_edges(res, min_k = min_k, fdr = fdr_cut, min_abs_r = min_abs_r, max_I2 = max_I2)
  cat(sprintf("  Consistent pairs identified: %d (%d concordant, %d discordant)\n",
              nrow(edges),
              sum(edges$relationship == "concordant"),
              sum(edges$relationship == "discordant")))

  # Annotate gene IDs with representative symbols if available
  if (nrow(edges) > 0 && !is.null(cluster_annotation) && "cluster_id" %in% colnames(cluster_annotation)) {
    symbol_col <- if ("Representative_Gene" %in% colnames(cluster_annotation)) {
      "Representative_Gene"
    } else if ("All_Genes_Locus_Tags" %in% colnames(cluster_annotation)) {
      "All_Genes_Locus_Tags"
    } else {
      NULL
    }
    if (!is.null(symbol_col)) {
      sym_map <- setNames(as.character(cluster_annotation[[symbol_col]]), as.character(cluster_annotation$cluster_id))
      edges$gene1_symbol <- sym_map[edges$gene1]
      edges$gene2_symbol <- sym_map[edges$gene2]
    }
  }

  out_edges_file <- file.path(out_dir, paste0("gene_pair_consistency_", mode, ".csv"))
  write.csv(edges, out_edges_file, row.names = FALSE)
  cat("  Saved gene pair consistency to:", out_edges_file, "\n")

  # --- 2. Concordant modules ---
  modules <- find_modules(res, min_k = min_k, max_I2 = max_I2,
                          cut_h = module_cut_h, min_size = module_min_sz)

  if (length(modules) > 0) {
    mod_df <- data.frame(cluster_id = names(modules), module = modules, stringsAsFactors = FALSE)
    if (!is.null(cluster_annotation) && "cluster_id" %in% colnames(cluster_annotation)) {
      mod_df <- merge(mod_df, cluster_annotation, by = "cluster_id", all.x = TRUE)
    }
    out_mod_file <- file.path(out_dir, paste0("concordant_modules_", mode, ".csv"))
    write.csv(mod_df, out_mod_file, row.names = FALSE)
    cat(sprintf("  Concordant modules found: %d (covering %d genes)\n", length(unique(mod_df$module)), nrow(mod_df)))
    cat("  Saved concordant modules to:", out_mod_file, "\n")

    # --- 3. Module-vs-Module Meta-analysis ---
    mm <- module_meta(prep_list, modules, max_I2 = max_I2)
    if (!is.null(mm)) {
      out_mm_file <- file.path(out_dir, paste0("module_pair_consistency_", mode, ".csv"))
      write.csv(mm, out_mm_file, row.names = FALSE)
      cat("  Saved module pair consistency to:", out_mm_file, "\n")

      disc <- subset(mm, relationship == "discordant" & consistent & padj < fdr_cut)
      if (nrow(disc) > 0) {
        cat("\n  Consistently DISCORDANT module pairs:\n")
        print(disc)
      } else {
        cat("  No module pairs met all consistency and significance criteria for discordance.\n")
      }
    }
  } else {
    cat("  No concordant modules formed with current thresholds.\n")
  }

  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  invisible(list(res = res, edges = edges, modules = modules))
}

# Standalone execution support: auto-run if sourced in a session where dataset_dds_list exists
if (sys.nframe() == 0 && exists("dataset_dds_list") && length(dataset_dds_list) > 0) {
  run_concordance_meta_analysis(
    dataset_dds_list   = dataset_dds_list,
    out_dir            = if (exists("out_dir")) out_dir else ".",
    cluster_annotation = if (exists("cluster_annotation")) cluster_annotation else NULL,
    mode               = if (exists("mode")) mode else "with_group",
    min_k              = if (exists("min_k")) min_k else 3,
    fdr_cut            = if (exists("fdr_cut")) fdr_cut else 0.05,
    min_abs_r          = if (exists("min_abs_r")) min_abs_r else 0.5,
    max_I2             = if (exists("max_I2")) max_I2 else 50,
    module_cut_h       = if (exists("module_cut_h")) module_cut_h else 0.6,
    module_min_sz      = if (exists("module_min_sz")) module_min_sz else 5
  )
}
