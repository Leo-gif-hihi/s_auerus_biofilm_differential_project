# Deep dive: trace the exact failure in wildmeta with metafor 5.2.1
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

# Step through estimate_null.rma manually
cat("=== Manual walkthrough of estimate_null.rma ===\n\n")

C_mat <- constrain_zero(1)(coef(m))
cat("C_mat:\n")
print(C_mat)

# 1. find_env
eval_env <- wildmeta:::find_env(m)
cat("\neval_env found:", !is.null(eval_env), "\n")

# 2. Check model$X
cat("\nclass(m$X):", class(m$X), "\n")
# In metafor 5.x, $ is overloaded. Use [[ instead
cat("class(m[['X']]):", class(m[["X"]]), "\n")
cat("m[['X']]:\n")
print(m[["X"]])

# 3. constrain_predictors
cat("\n=== constrain_predictors ===\n")
cat("Source:\n")
print(wildmeta:::constrain_predictors)

cat("\nCalling constrain_predictors(m$X, C_mat)...\n")
cat("  m$X class:", class(m$X), "\n")

# The issue: wildmeta uses m$X but metafor 5.x has $ overloading
# that makes m$X return the full model instead of the X matrix
tryCatch({
  res <- wildmeta:::constrain_predictors(m$X, C_mat)
  cat("  Result:\n")
  print(res)
}, error = function(e) {
  cat("  ERROR:", conditionMessage(e), "\n")
})

# Now try with m[["X"]]
cat("\nCalling constrain_predictors(m[['X']], C_mat)...\n")
tryCatch({
  res <- wildmeta:::constrain_predictors(m[["X"]], C_mat)
  cat("  Result:\n")
  print(res)
}, error = function(e) {
  cat("  ERROR:", conditionMessage(e), "\n")
})

# Check what m$X actually gives us
cat("\n=== What does m$X return? ===\n")
cat("is.matrix(m$X):", is.matrix(m$X), "\n")
cat("inherits(m$X, 'rma'):", inherits(m$X, "rma"), "\n")

# Check the $ method
cat("\n=== metafor $ method ===\n")
cat("Has $.rma.mv method:", existsMethod("$", "rma.mv") || !is.null(getS3method("$", "rma.mv", optional=TRUE)), "\n")
cat("Has $.rma method:", !is.null(getS3method("$", "rma", optional=TRUE)), "\n")

# Check if m$X is the same as m[["X"]]
cat("\nidentical(m$X, m[['X']]):", identical(m$X, m[["X"]]), "\n")
