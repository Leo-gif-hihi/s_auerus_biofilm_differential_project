library(limma)
set.seed(123)
y <- matrix(rnorm(100 * 8), nrow = 100, ncol = 8)
condition <- factor(rep(c("A", "B", "C", "D"), each = 2))
batch <- factor(rep(c("b1", "b2"), times = 4))
design <- model.matrix(~ 0 + condition + batch)
colnames(design)[1:4] <- c("A", "B", "C", "D")
fit <- lmFit(y, design)
contrast.matrix <- makeContrasts(B - A, C - A, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

unscaled_cov <- fit2$cov.coef
gene_vars <- fit2$s2.post
contrast_names <- colnames(fit2)
num_contrasts <- length(contrast_names)

cov_col_names <- c()
unscaled_cov_vals <- c()

for (i in 1:num_contrasts) {
  for (j in i:num_contrasts) {
    if (i == j) {
      cov_col_names <- c(cov_col_names, paste0("var_", contrast_names[i]))
    } else {
      cov_col_names <- c(cov_col_names, paste0("cov_", contrast_names[i], "_", contrast_names[j]))
    }
    unscaled_cov_vals <- c(unscaled_cov_vals, unscaled_cov[i, j])
  }
}

gene_cov_matrix <- outer(gene_vars, unscaled_cov_vals)
colnames(gene_cov_matrix) <- cov_col_names

print("unscaled_cov:")
print(unscaled_cov)
print("first row of gene_cov_matrix:")
print(gene_cov_matrix[1, ])

print("Manual calculation for gene 1:")
manual_cov <- unscaled_cov * gene_vars[1]
print(manual_cov)

# Check if fit2$cov.coef works
print("Does fit2$cov.coef work (partial matching)?")
print(!is.null(fit2$cov.coef))

# Compare the generated outer product covariance matrix with topTable (logFC / t)^2
tt <- topTable(fit2, coef = "B - A", number = Inf, sort.by = "none")
se_from_tstat <- tt$logFC[1] / tt$t[1]
var_from_tstat <- se_from_tstat^2

print("--- Comparison with logFC/t for Gene 1, contrast B-A ---")
print(paste("Variance from outer product matrix:", gene_cov_matrix[1, "var_B - A"]))
print(paste("Variance from (logFC/t)^2:", var_from_tstat))
print(paste("Are variances strictly equal?", all.equal(gene_cov_matrix[1, "var_B - A"], var_from_tstat)))

print("--- Simulation to verify covariance of estimated contrasts ---")
# Simulate data with a known variance (true_var) to verify that the empirical covariance
# of the estimated contrasts across many genes converges to true_var * fit2$cov.coef
set.seed(456)
n_genes <- 100000
true_var <- 1.5
y_sim <- matrix(rnorm(n_genes * 8, mean = 0, sd = sqrt(true_var)), nrow = n_genes, ncol = 8)

fit_sim <- lmFit(y_sim, design)
fit2_sim <- contrasts.fit(fit_sim, contrast.matrix)

# We use the unscaled covariance matrix from the fit (fit2_sim$cov.coef)
# The empirical covariance of the estimated contrasts should be equal to true_var * unscaled_cov
empirical_cov <- cov(fit2_sim$coefficients)
expected_cov <- true_var * fit2_sim$cov.coef

print("Empirical covariance of contrast estimates (across 100k simulated genes):")
print(empirical_cov)

print("Expected covariance (true_var * fit2$cov.coef):")
print(expected_cov)

print("Are empirical and expected approximately equal (tolerance 0.05)?")
print(all.equal(empirical_cov, expected_cov, tolerance = 0.05, check.attributes = FALSE))
