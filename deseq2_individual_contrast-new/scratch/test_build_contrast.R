library(DESeq2)
dds <- makeExampleDESeqDataSet(n=10, m=12)
dds$condition <- factor(rep(c("A", "B", "C"), each=4))
dds <- DESeq(dds, betaPrior=FALSE)
res_names <- resultsNames(dds)

build_contrast <- function(factor_name, num_level, den_level, res_names) {
  factor_terms <- grep(paste0("^", factor_name, "_"), res_names, value = TRUE)
  if (length(factor_terms) > 0) {
    base_level <- sub(".*_vs_", "", factor_terms[1])
    if (num_level == base_level) {
      c_num <- rep(0, length(res_names))
    } else {
      num_name <- paste0(factor_name, "_", num_level, "_vs_", base_level)
      c_num <- as.numeric(res_names == num_name)
    }
    if (den_level == base_level) {
      c_den <- rep(0, length(res_names))
    } else {
      den_name <- paste0(factor_name, "_", den_level, "_vs_", base_level)
      c_den <- as.numeric(res_names == den_name)
    }
    return(c_num - c_den)
  } else {
    stop("Could not find terms for factor")
  }
}
print(res_names)
print(build_contrast("condition", "B", "A", res_names))
print(build_contrast("condition", "C", "A", res_names))
print(build_contrast("condition", "C", "B", res_names))
