# Trace the math inside constrain_predictors
suppressPackageStartupMessages({
  library(metafor)
  library(clubSandwich)
  library(wildmeta)
})

gse <- c(rep("A", 3), rep("B", 2), rep("C", 3))
comp <- paste0("c", 1:8)
set.seed(42)
yi_vals <- rnorm(8, mean = 1, sd = 0.5)
vi_vals <- runif(8, 0.1, 0.5)

dat <- data.frame(yi = yi_vals, vi = vi_vals, gse = gse, comp = comp)
m <- rma.mv(yi = yi, V = vi, random = ~ 1 | gse / comp, data = dat, method = "REML")

C_mat <- constrain_zero(1)(coef(m))
Xmat <- m[["X"]]

cat("Xmat dimensions:", dim(Xmat), "\n")
cat("C_mat dimensions:", dim(C_mat), "\n")

q <- nrow(C_mat)
p <- ncol(C_mat)
cat("q (constraints):", q, "\n")
cat("p (coefficients):", p, "\n")
cat("p - q:", p - q, "\n")

# This is the issue: p=1, q=1, so p-q=0
# We're testing the intercept=0 on an intercept-only model
# That means Cnull_reduced has 0 columns -> Xnull has 0 columns -> can't multiply

# The constrain_zero(1) approach constrains ALL coefficients to zero.
# For an intercept-only model, the null model has no predictors at all!
# This causes the matrix multiplication to fail.

Cnull <- diag(nrow = p) - t(C_mat) %*% chol2inv(chol(tcrossprod(C_mat))) %*% C_mat
cat("\nCnull:\n")
print(Cnull)

svd_res <- svd(Cnull, nu = p - q, nv = p - q)
cat("\nSVD result:\n")
str(svd_res)

cat("\nCnull_reduced:\n")
print(svd_res$v)
cat("class:", class(svd_res$v), "\n")
cat("dim:", dim(svd_res$v), "\n")

cat("\n\nAttempting Xmat %*% Cnull_reduced:\n")
tryCatch({
  res <- Xmat %*% svd_res$v
  print(res)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat("This confirms: svd$v returns a 1x0 matrix when p-q=0\n")
  cat("Matrix multiplication 8x1 %*% 1x0 -> 'requires numeric/complex'\n")
})

cat("\n=== ROOT CAUSE ===\n")
cat("The issue is NOT a metafor version incompatibility.\n")
cat("It's a fundamental limitation: Wald_test_cwb with constrain_zero(1)\n")
cat("on an intercept-only model tries to create a null model with 0 predictors.\n")
cat("This is impossible because constrain_predictors does SVD on a rank-0\n")
cat("null space, producing a 0-column Cnull_reduced matrix.\n\n")
cat("The correct approach for an intercept-only model is NOT constrain_zero.\n")
cat("We need to use the Wald test differently.\n\n")

# Try clubSandwich::Wald_test directly (this should work fine)
cat("=== clubSandwich::Wald_test (direct, no bootstrapping) ===\n")
full_vcov <- vcovCR(m, type = "CR2", cluster = gse)
wt <- Wald_test(m, constraints = constrain_zero(1), vcov = full_vcov, test = "HTZ")
cat("Wald_test result:\n")
print(wt)

# The wildmeta CWB must use a different approach for intercept-only models
# It needs a "null model" to get residuals from, but for intercept-only
# the null model is yi = 0 (no predictors at all)

cat("\n=== Try: Fit model with a dummy moderator, then constrain ===\n")
dat$dummy <- rnorm(8)
m_mod <- rma.mv(yi = yi, V = vi, mods = ~ dummy, random = ~ 1 | gse / comp,
                data = dat, method = "REML")
cat("Model with dummy moderator fitted OK\n")
cat("Coefficients:", names(coef(m_mod)), "\n")

# Constrain intercept to zero (now there are 2 coefficients)
cwb_mod <- tryCatch({
  Wald_test_cwb(m_mod, constrain_zero(1), R = 49, adjust = "CR2", seed = 1)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(cwb_mod)) {
  cat("SUCCESS with moderator model!\n")
  print(cwb_mod)
}

cat("\n=== Try: constrain_zero for coefficient 2 only (dummy) ===\n")
cwb_mod2 <- tryCatch({
  Wald_test_cwb(m_mod, constrain_zero(2), R = 49, adjust = "CR2", seed = 1)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(cwb_mod2)) {
  cat("SUCCESS constraining coefficient 2!\n")
  print(cwb_mod2)
}
