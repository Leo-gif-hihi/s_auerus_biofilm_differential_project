library(DESeq2)
dds <- makeExampleDESeqDataSet(n=3000, m=12)
dds$condition <- factor(rep(c("A", "B", "C"), each=4))
dds <- DESeq(dds, betaPrior=FALSE)

t0 <- Sys.time()
cov_list <- lapply(1:nrow(dds), function(i) {
  DESeq2:::covarianceMatrix(dds, rowNumber=i)
})
t1 <- Sys.time()
print(t1 - t0)

resNames <- resultsNames(dds)
print(resNames)

# Get contrast numeric vector
c_B <- DESeq2:::makeContrast(c("condition", "B", "A"), resNames)
c_C <- DESeq2:::makeContrast(c("condition", "C", "A"), resNames)
print(c_B)
print(c_C)

# Covariance
cov_mat <- cov_list[[1]]
var_B <- t(c_B) %*% cov_mat %*% c_B
var_C <- t(c_C) %*% cov_mat %*% c_C
cov_BC <- t(c_B) %*% cov_mat %*% c_C

print(c(var_B, var_C, cov_BC))
