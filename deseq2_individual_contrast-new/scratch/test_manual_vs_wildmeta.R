# ============================================================================
# COMPARING MANUAL CWB VS WILDMETA PACKAGE
# ============================================================================
# Since wildmeta fails on intercept-only models (which we use), we will
# test a model with a single moderator (where wildmeta works). 
# We'll run BOTH wildmeta and our manual CWB approach on this model to prove
# that our manual approach produces the same results as the package.

suppressPackageStartupMessages({
  library(metafor)
  library(clubSandwich)
  library(wildmeta)
})

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("TEST: Manual CWB vs wildmeta::Wald_test_cwb\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# 1. Create a synthetic dataset
set.seed(42)
gse_full <- c(rep("GSE001", 3), rep("GSE002", 2), rep("GSE003", 4), rep("GSE004", 3), rep("GSE005", 2))
comp_full <- paste0("comp_", seq_along(gse_full))
n_comp <- length(gse_full)

yi_val <- rnorm(n_comp, mean = 0.5, sd = 0.8)
vi_val <- runif(n_comp, 0.15, 0.6)
mod1 <- rnorm(n_comp) # Dummy moderator so wildmeta works

dat <- data.frame(yi = yi_val, vi = vi_val, gse = gse_full, comp = comp_full, mod1 = mod1)
V_mat <- vcalc(vi = vi_val, cluster = gse_full, obs = comp_full, rho = 0.5)

# 2. Fit the Full Model
cat("Fitting full model with moderator...\n")
full_mod <- rma.mv(yi = yi, V = V_mat, mods = ~ mod1, random = ~ 1 | gse / comp, data = dat, method = "REML")

# 3. Wildmeta CWB (Testing moderator = 0)
cat("\nRunning wildmeta::Wald_test_cwb (adjust = 'CR0', Webb6 distribution)...\n")
set.seed(123)
res_wild <- Wald_test_cwb(
  full_model = full_mod,
  constraints = constrain_zero(2), # Test moderator (coefficient 2)
  R = 399,
  adjust = "CR0", # Matches our manual implementation (raw residuals)
  type = "CR2",
  auxiliary_dist = "Webb six",
  seed = 123
)

wild_F <- attr(res_wild, "original")
wild_p <- res_wild$p_val
cat(sprintf("  wildmeta F-statistic: %.4f\n", wild_F))
cat(sprintf("  wildmeta p-value:     %.4f\n", wild_p))

# Extract the run_manual_cwb function from the main pipeline script
# We parse the file and evaluate only the expression that defines the function
cat("Extracting run_manual_cwb from DESeq2_individual_analysis-new.R...\n")
pipeline_exprs <- parse("DESeq2_individual_analysis-new.R")
for (e in pipeline_exprs) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("run_manual_cwb"))) {
    eval(e)
    break
  }
}

if (!exists("run_manual_cwb")) {
  stop("Failed to extract run_manual_cwb function from pipeline script.")
} else {
  cat("Successfully imported run_manual_cwb.\n")
}

# 4. Manual CWB (Testing moderator = 0)
cat("\nRunning Manual CWB approach...\n")

# Get original F-statistic using clubSandwich
orig_robust <- robust(full_mod, cluster = gse_full, clubSandwich = TRUE)
orig_F_manual <- as.numeric(orig_robust$beta[2])^2 / as.numeric(orig_robust$se[2])^2

# Fit Null model (without moderator) to get residuals
null_mod <- rma.mv(yi = yi, V = V_mat, random = ~ 1 | gse / comp, data = dat, method = "REML")
preds <- fitted(null_mod)
resids <- yi_val - preds

# Define the refit function
refit_fn <- function(y_boot) {
  dat_boot <- dat
  dat_boot$yi <- y_boot
  boot_mod <- tryCatch(
    rma.mv(yi = yi, V = V_mat, mods = ~ mod1, random = ~ 1 | gse / comp, data = dat_boot, method = "REML"),
    error = function(e) NULL
  )
  if (!is.null(boot_mod)) {
    rob_boot <- tryCatch(
      robust(boot_mod, cluster = gse_full, clubSandwich = TRUE),
      error = function(e) NULL
    )
    if (!is.null(rob_boot)) {
      return(as.numeric(rob_boot$beta[2])^2 / as.numeric(rob_boot$se[2])^2)
    }
  }
  return(NA)
}

cwb_res <- run_manual_cwb(
  orig_yi = yi_val,
  orig_F = orig_F_manual,
  cluster_var = gse_full,
  refit_fn = refit_fn,
  R = 399,
  null_preds = preds,
  null_resids = resids,
  seed = 123
)

manual_p <- cwb_res$p_val

cat(sprintf("  Manual F-statistic: %.4f\n", orig_F_manual))
cat(sprintf("  Manual p-value:     %.4f\n", manual_p))

cat("\n======================================================================\n")
cat("CONCLUSION\n")
cat("======================================================================\n")
if (abs(wild_p - manual_p) < 0.05) {
  cat("SUCCESS: The manual approach and the wildmeta package produce\n")
  cat("statistically equivalent results!\n")
} else {
  cat("DIFFERENCE DETECTED: The results do not match.\n")
}
cat("======================================================================\n")
