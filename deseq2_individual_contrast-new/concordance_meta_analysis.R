# ============================================================================
# Consistently CONCORDANT / DISCORDANT gene groups across datasets
# (random-effects meta-analysis of gene-gene association)
# ============================================================================
#
# HOW TO RUN
#   Source this file AFTER the DESeq2 loop in DESeq2_individual_analysis-new.R,
#   i.e. when `dataset_dds_list` (one DESeqDataSet per GSE) and `out_dir`
#   exist in the session. (You can also `source()` it from the end of that script.)
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

# ---- settings --------------------------------------------------------------
mode          <- "with_group"  # "with_group" or "within_group"
min_k         <- 3             # min # datasets in which BOTH genes are measured
fdr_cut       <- 0.05
min_abs_r     <- 0.5           # min |pooled r| to call an edge
max_I2        <- 50            # max heterogeneity (%) to call "consistent"
module_cut_h  <- 0.6           # hclust cut height on (1 - r); 0.6 ~ avg r >= 0.4
module_min_sz <- 5

if (!exists("out_dir")) out_dir <- "."

fisher_z <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))

# ---- 1. per-dataset expression (VST, batch removed) ------------------------
prep_expr <- function(dds, mode = c("with_group", "within_group")) {
  mode  <- match.arg(mode)
  cd    <- as.data.frame(colData(dds))
  vsd   <- tryCatch(
    assay(vst(dds, blind = FALSE)),
    error = function(e) assay(varianceStabilizingTransformation(dds, blind = FALSE))
  )
  batch <- droplevels(factor(cd$Batch))
  group <- droplevels(factor(cd$Group))

  if (mode == "with_group") {
    # remove batch, protect the group effect
    x <- if (nlevels(batch) > 1) {
      removeBatchEffect(vsd, batch = batch, design = model.matrix(~group))
    } else vsd
    p <- nlevels(batch) - 1
  } else {
    # remove batch AND group (residuals)
    des <- if (nlevels(batch) > 1) model.matrix(~ batch + group) else model.matrix(~group)
    x   <- residuals(lmFit(vsd, des), vsd)
    p   <- ncol(des) - 1
  }
  # df = n - p - 3 is the inverse of Var(Fisher z) for a partial correlation
  list(x = x, df = ncol(x) - p - 3)
}

# ---- 2. all-pairs random-effects meta-analysis (vectorised, 2 passes) ------
coexpr_meta <- function(prep_list, min_k = 3, use_hk = TRUE) {
  prep_list <- prep_list[sapply(prep_list, function(p) p$df >= 1)]
  genes <- sort(unique(unlist(lapply(prep_list, function(p) rownames(p$x)))))
  G     <- length(genes)
  cat("Genes:", G, "| datasets used:", length(prep_list), "\n")
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
    K  <- K  + ok
    W  <- W  + w * ok
    W2 <- W2 + w^2 * ok
    WZ <- WZ + w * Z
    WZ2 <- WZ2 + w * Z^2
    NP <- NP + (Z > 0)
    NN <- NN + (Z < 0)
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
    ws <- ok / (1 / p$df + tau2)
    Ws   <- Ws   + ws
    WsZ  <- WsZ  + ws * Z
    WsZ2 <- WsZ2 + ws * Z^2
  }
  zhat <- WsZ / Ws
  if (use_hk) {
    qhk <- (WsZ2 - WsZ^2 / Ws) / (K - 1)
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
  sel <- upper.tri(res$r) & res$K >= min_k & !is.na(res$p)
  idx <- which(sel)
  padj <- p.adjust(res$p[idx], method = "BH")        # BH over ALL tested pairs
  same_sign <- res$NP[idx] == 0 | res$NN[idx] == 0   # same direction in every dataset
  pass <- padj < fdr & abs(res$r[idx]) >= min_abs_r & res$I2[idx] <= max_I2 & same_sign
  padj <- padj[pass]; idx <- idx[pass]
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
  R <- res$r
  bad <- is.na(R) | res$K < min_k | res$I2 > max_I2 | (res$NP > 0 & res$NN > 0)
  R[bad] <- 0
  diag(R) <- 1
  keep <- rowSums(abs(R) >= 0.5) > 1
  R <- R[keep, keep, drop = FALSE]
  hc <- hclust(as.dist(1 - R), method = "average")
  cl <- cutree(hc, h = cut_h)
  tab <- table(cl)
  good <- as.integer(names(tab)[tab >= min_size])
  cl[!cl %in% good] <- NA
  lab <- setNames(paste0("M", match(cl, good)), names(cl))
  lab[is.na(cl)] <- NA
  lab
}

# ---- 5. module-vs-module random-effects test (few tests, less noise) -------
module_meta <- function(prep_list, modules) {
  mods <- sort(unique(na.omit(modules)))
  if (length(mods) < 2) return(NULL)
  per_ds <- lapply(prep_list, function(p) {
    if (p$df < 1) return(NULL)
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
    pi <- predict(fit)
    out[[length(out) + 1]] <- data.frame(
      module1 = mods[i], module2 = mods[j], k_datasets = sum(ok),
      pooled_r = tanh(as.numeric(fit$beta)),
      ci_lb = tanh(fit$ci.lb), ci_ub = tanh(fit$ci.ub),
      pi_lb = tanh(pi$pi.lb), pi_ub = tanh(pi$pi.ub),
      pval = fit$pval, I2 = fit$I2, tau2 = fit$tau2,
      n_datasets_pos = sum(z[ok] > 0), n_datasets_neg = sum(z[ok] < 0)
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
# RUN
# ============================================================================
cat("\n=== GENE-PAIR CONCORDANCE META-ANALYSIS (mode:", mode, ") ===\n")
prep_list <- lapply(dataset_dds_list, prep_expr, mode = mode)
print(sapply(prep_list, function(p) p$df))  # effective n - 3 per dataset; need >= 1

res <- coexpr_meta(prep_list, min_k = min_k)

edges <- call_edges(res, min_k = min_k, fdr = fdr_cut, min_abs_r = min_abs_r, max_I2 = max_I2)
cat("Consistent pairs:", nrow(edges), "(",
    sum(edges$relationship == "concordant"), "concordant,",
    sum(edges$relationship == "discordant"), "discordant )\n")
write.csv(edges, file.path(out_dir, paste0("gene_pair_consistency_", mode, ".csv")), row.names = FALSE)

modules <- find_modules(res, min_k = min_k, max_I2 = max_I2,
                        cut_h = module_cut_h, min_size = module_min_sz)
mod_df <- data.frame(cluster_id = names(modules), module = modules, stringsAsFactors = FALSE)
mod_df <- mod_df[!is.na(mod_df$module), ]
if (exists("cluster_annotation")) {
  mod_df$cluster_id_key <- mod_df$cluster_id
  mod_df <- merge(mod_df, cluster_annotation, by = "cluster_id", all.x = TRUE)
  mod_df$cluster_id_key <- NULL
}
write.csv(mod_df, file.path(out_dir, paste0("concordant_modules_", mode, ".csv")), row.names = FALSE)
cat("Modules found:", length(unique(mod_df$module)), "\n")

mm <- module_meta(prep_list, modules)
if (!is.null(mm)) {
  write.csv(mm, file.path(out_dir, paste0("module_pair_consistency_", mode, ".csv")), row.names = FALSE)
  cat("\nModule pairs that are consistently DISCORDANT:\n")
  print(subset(mm, relationship == "discordant" & consistent & padj < 0.05))
}

# CAVEAT: modules are learned from the same datasets they are then tested in,
# so module-level p-values are optimistic. For a clean test, learn modules on
# K-1 datasets and run module_meta() on the held-out one (leave-one-dataset-out).
