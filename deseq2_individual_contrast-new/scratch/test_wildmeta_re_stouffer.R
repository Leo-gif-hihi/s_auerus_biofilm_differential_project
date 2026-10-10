# ============================================================================
# TEST: Does wildmeta::Wald_test_cwb work with our RE Stouffer model?
# ============================================================================
#
# This script creates a small synthetic dataset that mimics the exact structure
# of the RE Stouffer model in DESeq2_individual_analysis-new.R:
#
#   metafor::rma.mv(
#     yi = z_i_twosided,        # Two-sided Z-scores from p-values
#     V  = Cor_mat_stouffer,     # Empirical correlation matrix
#     W  = neff_vec,             # Effective sample size weights
#     random = ~ 1 | GSE / comparison,  # 3-level RE structure
#     method = "REML"
#   )
#
# We then apply metafor::robust() + wildmeta::Wald_test_cwb() to test
# whether CWB is compatible with this model specification.
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("TEST: wildmeta CWB with RE Stouffer model (small synthetic data)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# --- 1. Load required packages ---
suppressPackageStartupMessages({
  library(metafor)
  library(clubSandwich)
  library(wildmeta)
})

cat("Packages loaded:\n")
cat("  metafor:      ", as.character(packageVersion("metafor")), "\n")
cat("  clubSandwich: ", as.character(packageVersion("clubSandwich")), "\n")
cat("  wildmeta:     ", as.character(packageVersion("wildmeta")), "\n\n")

# --- 2. Create synthetic data mimicking our pipeline ---
# Simulate 5 GSE studies, some with multiple comparisons (14 comparisons total)
set.seed(42)

# Study structure: GSE1 has 3 comparisons, GSE2 has 2, GSE3 has 4, GSE4 has 3, GSE5 has 2
gse_ids <- c(
  rep("GSE001", 3),
  rep("GSE002", 2),
  rep("GSE003", 4),
  rep("GSE004", 3),
  rep("GSE005", 2)
)
n_comparisons <- length(gse_ids)

comparison_ids <- paste0("comp_", seq_len(n_comparisons))

# Effective sample sizes (typical range from your data)
neff_vec <- c(
  2.5, 3.0, 2.0,      # GSE001
  4.0, 3.5,            # GSE002
  2.0, 2.5, 3.0, 2.0,  # GSE003
  3.5, 4.0, 3.0,       # GSE004
  2.5, 3.0              # GSE005
)
names(neff_vec) <- comparison_ids

cat("Study structure:\n")
for (g in unique(gse_ids)) {
  idx <- which(gse_ids == g)
  cat(sprintf("  %s: %d comparisons (Neff: %s)\n",
              g, length(idx), paste(neff_vec[idx], collapse = ", ")))
}

# --- 3. Simulate p-values and convert to Z-scores (like our pipeline) ---
# Simulate a gene with a real signal (p-values mostly small)
set.seed(123)
p_vals <- c(
  0.01, 0.03, 0.15,     # GSE001
  0.005, 0.02,           # GSE002
  0.08, 0.04, 0.12, 0.2, # GSE003
  0.001, 0.01, 0.05,     # GSE004
  0.03, 0.06              # GSE005
)
names(p_vals) <- comparison_ids

# Simulate LFC signs (mostly positive = consistent direction)
lfc_signs <- c(1, 1, -1, 1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1)

# Two-sided Z-scores (what our pipeline uses as yi)
z_twosided <- qnorm(1 - p_vals)

# Direction-aware Z-scores (for computing correlation matrix)
z_directed <- lfc_signs * qnorm(1 - p_vals / 2)

cat("\n\nTwo-sided Z-scores (yi):\n")
print(round(z_twosided, 3))

# --- 4. Build empirical correlation matrix (like our pipeline) ---
# In the real pipeline, this comes from cor(Z_matrix, use="pairwise.complete.obs")
# Here we construct a block-diagonal correlation matrix with within-GSE correlations

Cor_mat <- diag(1, n_comparisons)
rownames(Cor_mat) <- comparison_ids
colnames(Cor_mat) <- comparison_ids

# Add within-GSE correlations (rho = 0.5 as used in vcalc)
for (g in unique(gse_ids)) {
  idx <- which(gse_ids == g)
  if (length(idx) > 1) {
    for (r in idx) {
      for (cc in idx) {
        if (r != cc) {
          Cor_mat[r, cc] <- 0.5
        }
      }
    }
  }
}

cat("\nCorrelation matrix (V):\n")
print(round(Cor_mat, 2))

# --- 5. Fit the RE Stouffer model (exact same specification as pipeline) ---
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("FITTING RE STOUFFER MODEL (rma.mv)\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

res_rma <- tryCatch({
  rma.mv(
    yi     = z_twosided,
    V      = Cor_mat,
    W      = neff_vec,
    random = ~ 1 | gse_ids / comparison_ids,
    method = "REML"
  )
}, error = function(e) {
  cat("ERROR fitting rma.mv:", e$message, "\n")
  NULL
})

if (!is.null(res_rma)) {
  cat("rma.mv fitted successfully!\n")
  cat("  Beta (pooled Z):  ", round(as.numeric(res_rma$beta), 4), "\n")
  cat("  SE:               ", round(as.numeric(res_rma$se), 4), "\n")
  cat("  p-value:          ", format(res_rma$pval, digits = 4), "\n")
  cat("  sigma2 (Level 3): ", round(res_rma$sigma2[1], 4), "\n")
  cat("  sigma2 (Level 2): ", round(res_rma$sigma2[2], 4), "\n")
} else {
  cat("rma.mv FAILED. Cannot proceed.\n")
  quit(save = "no", status = 1)
}

# --- 6. Apply RVE (robust variance estimation) ---
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("APPLYING ROBUST VARIANCE ESTIMATION (clubSandwich)\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

robust_res <- tryCatch({
  robust(res_rma, cluster = gse_ids, clubSandwich = TRUE)
}, error = function(e) {
  cat("ERROR in robust():", e$message, "\n")
  NULL
})

if (!is.null(robust_res)) {
  cat("robust() applied successfully!\n")
  cat("  Robust Beta:  ", round(as.numeric(robust_res$beta), 4), "\n")
  cat("  Robust SE:    ", round(as.numeric(robust_res$se), 4), "\n")
  cat("  Robust p-val: ", format(robust_res$pval, digits = 4), "\n")
  cat("  Robust ddf:   ", round(robust_res$ddf, 2), "\n")
} else {
  cat("robust() FAILED. Cannot proceed.\n")
  quit(save = "no", status = 1)
}

# --- 7. TEST: wildmeta::Wald_test_cwb() ---
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TESTING wildmeta::Wald_test_cwb()\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Test with different adjustment methods and working models
adjustment_methods <- c("CR0", "CR1", "CR2")

for (adj in adjustment_methods) {
  cat(sprintf("  Testing: adjust = '%s' ... ", adj))

  cwb_result <- tryCatch({
    Wald_test_cwb(
      full_model   = res_rma,
      constraints  = constrain_zero(1),  # Test intercept = 0
      R            = 99,                 # Small R for speed
      adjust       = adj,
      type         = "CWB",
      seed         = 42
    )
  }, error = function(e) {
    list(error = e$message)
  })

  if ("error" %in% names(cwb_result)) {
    cat("FAILED:", cwb_result$error, "\n")
  } else {
    cat(sprintf("p_boot = %.4f (Fstat = %.3f)\n",
                cwb_result$p_val, cwb_result$Fstat))
  }
}

# --- 8. Detailed test with our preferred settings (CR2 + more reps) ---
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("DETAILED TEST: CR2 + R=399 (our target settings)\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

cwb_detailed <- tryCatch({
  Wald_test_cwb(
    full_model   = res_rma,
    constraints  = constrain_zero(1),
    R            = 399,
    adjust       = "CR2",
    type         = "CWB",
    seed         = 12345
  )
}, error = function(e) {
  list(error = e$message)
})

if ("error" %in% names(cwb_detailed)) {
  cat("FAILED:", cwb_detailed$error, "\n\n")
} else {
  cat("Wald_test_cwb() SUCCEEDED with CR2 + R=399!\n\n")
  cat("  Results:\n")
  print(cwb_detailed)
  cat("\n")
}

# --- 9. Edge case: Only 2 GSE studies (minimum for meta-analysis) ---
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("EDGE CASE: 2 GSE studies only\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

idx_2gse <- gse_ids %in% c("GSE001", "GSE002")
z_small <- z_twosided[idx_2gse]
gse_small <- gse_ids[idx_2gse]
comp_small <- comparison_ids[idx_2gse]
neff_small <- neff_vec[idx_2gse]
Cor_small <- Cor_mat[idx_2gse, idx_2gse]

res_rma_2gse <- tryCatch({
  rma.mv(
    yi     = z_small,
    V      = Cor_small,
    W      = neff_small,
    random = ~ 1 | gse_small / comp_small,
    method = "REML"
  )
}, error = function(e) {
  cat("  rma.mv with 2 GSE: ERROR:", e$message, "\n")
  NULL
})

if (!is.null(res_rma_2gse)) {
  cat("  rma.mv with 2 GSE: fitted\n")

  cwb_2gse <- tryCatch({
    Wald_test_cwb(
      full_model   = res_rma_2gse,
      constraints  = constrain_zero(1),
      R            = 99,
      adjust       = "CR2",
      type         = "CWB",
      seed         = 42
    )
  }, error = function(e) {
    list(error = e$message)
  })

  if ("error" %in% names(cwb_2gse)) {
    cat("  CWB with 2 GSE:   FAILED:", cwb_2gse$error, "\n")
  } else {
    cat(sprintf("  CWB with 2 GSE:   p_boot = %.4f\n", cwb_2gse$p_val))
  }
}

# --- 10. Edge case: 3 GSE studies ---
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("EDGE CASE: 3 GSE studies\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

idx_3gse <- gse_ids %in% c("GSE001", "GSE002", "GSE003")
z_3 <- z_twosided[idx_3gse]
gse_3 <- gse_ids[idx_3gse]
comp_3 <- comparison_ids[idx_3gse]
neff_3 <- neff_vec[idx_3gse]
Cor_3 <- Cor_mat[idx_3gse, idx_3gse]

res_rma_3gse <- tryCatch({
  rma.mv(
    yi     = z_3,
    V      = Cor_3,
    W      = neff_3,
    random = ~ 1 | gse_3 / comp_3,
    method = "REML"
  )
}, error = function(e) {
  cat("  rma.mv with 3 GSE: ERROR:", e$message, "\n")
  NULL
})

if (!is.null(res_rma_3gse)) {
  cat("  rma.mv with 3 GSE: fitted\n")

  cwb_3gse <- tryCatch({
    Wald_test_cwb(
      full_model   = res_rma_3gse,
      constraints  = constrain_zero(1),
      R            = 99,
      adjust       = "CR2",
      type         = "CWB",
      seed         = 42
    )
  }, error = function(e) {
    list(error = e$message)
  })

  if ("error" %in% names(cwb_3gse)) {
    cat("  CWB with 3 GSE:   FAILED:", cwb_3gse$error, "\n")
  } else {
    cat(sprintf("  CWB with 3 GSE:   p_boot = %.4f\n", cwb_3gse$p_val))
  }
}

# --- 11. Edge case: Gene with no signal (null gene) ---
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("EDGE CASE: Null gene (no real signal)\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

set.seed(999)
p_null <- runif(n_comparisons, 0.3, 0.95)
z_null <- qnorm(1 - p_null)

res_rma_null <- tryCatch({
  rma.mv(
    yi     = z_null,
    V      = Cor_mat,
    W      = neff_vec,
    random = ~ 1 | gse_ids / comparison_ids,
    method = "REML"
  )
}, error = function(e) {
  cat("  rma.mv null gene: ERROR:", e$message, "\n")
  NULL
})

if (!is.null(res_rma_null)) {
  cat(sprintf("  Null gene rma.mv: beta = %.3f, p = %.4f\n",
              as.numeric(res_rma_null$beta), res_rma_null$pval))

  cwb_null <- tryCatch({
    Wald_test_cwb(
      full_model   = res_rma_null,
      constraints  = constrain_zero(1),
      R            = 99,
      adjust       = "CR2",
      type         = "CWB",
      seed         = 42
    )
  }, error = function(e) {
    list(error = e$message)
  })

  if ("error" %in% names(cwb_null)) {
    cat("  CWB null gene:    FAILED:", cwb_null$error, "\n")
  } else {
    cat(sprintf("  CWB null gene:    p_boot = %.4f (should be large)\n", cwb_null$p_val))
  }
}

# --- 12. Compare: RVE p-value vs wildmeta CWB p-value ---
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("COMPARISON: RVE p-value vs wildmeta CWB p-value\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Note: Our current pipeline manually implements CWB for the 3LRE model.\n")
cat("wildmeta automates this with additional corrections (e.g., CR2 adjustment).\n\n")

if (!("error" %in% names(cwb_detailed))) {
  cat(sprintf("  wildmeta CWB p-value: %.4f\n", cwb_detailed$p_val))
  cat(sprintf("  Robust RVE p-value:   %.4f\n", as.numeric(robust_res$pval)))
  cat("\n")
}

# --- 13. Test with the EXACT same model spec as the pipeline's 3LRE ---
cat("\n", paste(rep("-", 60), collapse = ""), "\n")
cat("BONUS: Testing CWB with 3LRE model (LFC-based, same as pipeline)\n")
cat(paste(rep("-", 60), collapse = ""), "\n\n")

# Simulate LFC values and SE (like the 3LRE model)
set.seed(77)
yi_lfc <- rnorm(n_comparisons, mean = 0.5, sd = 0.8)
sei_lfc <- runif(n_comparisons, 0.15, 0.6)
vi_lfc <- sei_lfc^2

V_mat_3lre <- metafor::vcalc(vi = vi_lfc, cluster = gse_ids,
                              obs = comparison_ids, rho = 0.5)

res_3lre <- tryCatch({
  rma.mv(
    yi     = yi_lfc,
    V      = V_mat_3lre,
    random = ~ 1 | gse_ids / comparison_ids,
    method = "REML"
  )
}, error = function(e) {
  cat("  3LRE rma.mv: ERROR:", e$message, "\n")
  NULL
})

if (!is.null(res_3lre)) {
  cat(sprintf("  3LRE rma.mv: beta = %.3f, SE = %.3f, p = %.4f\n",
              as.numeric(res_3lre$beta), as.numeric(res_3lre$se), res_3lre$pval))

  cwb_3lre <- tryCatch({
    Wald_test_cwb(
      full_model   = res_3lre,
      constraints  = constrain_zero(1),
      R            = 399,
      adjust       = "CR2",
      type         = "CWB",
      seed         = 42
    )
  }, error = function(e) {
    list(error = e$message)
  })

  if ("error" %in% names(cwb_3lre)) {
    cat("  CWB 3LRE:    FAILED:", cwb_3lre$error, "\n")
  } else {
    cat(sprintf("  CWB 3LRE:    p_boot = %.4f (Fstat = %.3f)\n",
                cwb_3lre$p_val, cwb_3lre$Fstat))
  }
}

# --- Summary ---
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat("If all tests above show p_boot values, then wildmeta::Wald_test_cwb()\n")
cat("is fully compatible with both the RE Stouffer and 3LRE model specs.\n")
cat("You can replace the manual CWB loop with wildmeta for:\n")
cat("  - Proper CR2 small-sample correction\n")
cat("  - Webb six-point distribution (same as your manual implementation)\n")
cat("  - Cleaner code and fewer potential bugs\n\n")
