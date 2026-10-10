library(DESeq2)

# Create a dummy dataset
dds <- makeExampleDESeqDataSet(n=100, m=6)
dds <- DESeq(dds, betaPrior=FALSE)

# Try to extract the covariance matrix for the first gene (row 1)
cov_mat <- DESeq2:::covarianceMatrix(dds, rowNumber=1)

print("Variance-covariance matrix for gene 1:")
print(cov_mat)

# Extract results for the same dummy dataset
res <- results(dds)

# The 'results' function gives the Log2 Fold Change Standard Error (lfcSE) 
# for the last coefficient by default (in this case, 'condition_B_vs_A' or 'conditionB')
# Let's compare the square root of the variance for 'conditionB' with the lfcSE from results()

lfcSE_from_results <- res$lfcSE[1]
se_from_cov_mat <- sqrt(cov_mat["conditionB", "conditionB"])

print("---")
print(paste("lfcSE from results():", lfcSE_from_results))
print(paste("Standard Error from sqrt of variance:", se_from_cov_mat))
print(paste("Are they exactly equal?", all.equal(lfcSE_from_results, se_from_cov_mat)))

# Also test with betaPrior=TRUE (default for some older versions, but false is default since 1.16)
tryCatch({
  dds_prior <- DESeq(dds, betaPrior=TRUE)
  cov_mat_prior <- DESeq2:::covarianceMatrix(dds_prior, rowNumber=1)
  print("Variance-covariance matrix for gene 1 (betaPrior=TRUE):")
  print(cov_mat_prior)
}, error = function(e) {
  print("Error with betaPrior=TRUE:")
  print(e$message)
})
