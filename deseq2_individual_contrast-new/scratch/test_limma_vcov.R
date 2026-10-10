library(limma)

# Create dummy data (100 genes, 6 samples)
set.seed(123)
y <- matrix(rnorm(100*6), nrow=100, ncol=6)
condition <- factor(rep(c("A", "B"), each=3))
design <- model.matrix(~ condition)

# Run limma pipeline
fit <- lmFit(y, design)
fit <- eBayes(fit)

# Extract results for the first gene (row 1)
# 1. Unscaled covariance matrix of the coefficients (same for all genes if no missing data/weights)
unscaled_cov_mat <- fit$cov.coefficients

# 2. Gene-specific posterior residual variance for gene 1
s2_post_gene1 <- fit$s2.post[1]

# 3. Scaled Variance-Covariance Matrix for gene 1
cov_mat_gene1 <- unscaled_cov_mat * s2_post_gene1

print("Variance-covariance matrix for gene 1 in limma:")
print(cov_mat_gene1)

# Compare with standard errors
# limma calculates standard error for the coefficients as:
# stdev.unscaled * sqrt(s2.post)
se_from_limma <- fit$stdev.unscaled[1, "conditionB"] * sqrt(s2_post_gene1)

# Standard error from our calculated covariance matrix
se_from_cov_mat <- sqrt(cov_mat_gene1["conditionB", "conditionB"])

# Standard error and variance calculated from logFC and t-statistic
tt <- topTable(fit, coef="conditionB", number=Inf, sort.by="none")
se_from_tstat <- tt$logFC[1] / tt$t[1]
var_from_tstat <- se_from_tstat^2

print("---")
print(paste("Standard Error from limma internals:", se_from_limma))
print(paste("Standard Error from sqrt of variance:", se_from_cov_mat))
print(paste("Standard Error from logFC/t:", se_from_tstat))
print(paste("Are SEs exactly equal (internal vs cov mat)?", all.equal(se_from_limma, se_from_cov_mat)))
print(paste("Are SEs exactly equal (logFC/t vs cov mat)?", all.equal(se_from_tstat, se_from_cov_mat)))
print("---")
print(paste("Variance from covariance matrix:", cov_mat_gene1["conditionB", "conditionB"]))
print(paste("Variance from (logFC/t)^2:", var_from_tstat))
print(paste("Are variances equal?", all.equal(var_from_tstat, cov_mat_gene1["conditionB", "conditionB"])))
