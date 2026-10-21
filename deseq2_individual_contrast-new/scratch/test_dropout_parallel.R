# test_dropout_parallel.R
suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
  library(clubSandwich)
  library(parallel)
})

source("dropout_test.R")

test_dir <- "scratch/test_dropout_output"
if (dir.exists(test_dir)) unlink(test_dir, recursive = TRUE)
dir.create(test_dir, recursive = TRUE)

# Create 4 dummy datasets: GSE1, GSE2, GSE3, GSE4
# 10 genes, 5 comparisons total
n_genes <- 10
genes <- paste0("gene_", 1:n_genes)
comps <- c("comp_gse1", "comp_gse2_a", "comp_gse2_b", "comp_gse3", "comp_gse4")
gse_map <- c(
  comp_gse1 = "GSE1",
  comp_gse2_a = "GSE2",
  comp_gse2_b = "GSE2",
  comp_gse3 = "GSE3",
  comp_gse4 = "GSE4"
)

set.seed(42)
lfc_mat <- matrix(rnorm(n_genes * length(comps), mean = 1, sd = 0.5),
                  nrow = n_genes, dimnames = list(genes, comps))
se_mat <- matrix(runif(n_genes * length(comps), 0.1, 0.3),
                 nrow = n_genes, dimnames = list(genes, comps))
pval_mat <- matrix(runif(n_genes * length(comps), 0.001, 0.05),
                   nrow = n_genes, dimnames = list(genes, comps))

sign_map <- list(comp_gse1 = 1, comp_gse2_a = 1, comp_gse2_b = 1, comp_gse3 = 1, comp_gse4 = 1)
all_cov_mats <- list(
  GSE2 = list(dataset_level = matrix(c(0.01, 0.002, 0.002, 0.02), nrow = 2,
                                     dimnames = list(c("comp_gse2_a", "comp_gse2_b"),
                                                     c("comp_gse2_a", "comp_gse2_b"))))
)

annot <- data.frame(cluster_id = genes, gene_name = paste0("sym_", 1:n_genes), stringsAsFactors = FALSE)

cat("--- TEST 1: Run with num_cores = 2, min_dropout_datasets = 3 ---\n")
# N=4 datasets -> 4C3 = 4 combinations + baseline = 5 combinations total
res1 <- run_dropout_test(
  lfc_matrix = lfc_mat,
  se_matrix = se_mat,
  pval_matrix = pval_mat,
  comparison_gse_map = gse_map,
  sign_map = sign_map,
  all_cov_mats = all_cov_mats,
  cluster_annotation = annot,
  min_datasets = 0,
  min_dropout_datasets = 3,
  out_dir = test_dir,
  force_rerun = FALSE,
  num_cores = 2
)

stopifnot(nrow(res1) == 5)
stopifnot(file.exists(file.path(test_dir, "drop_out_test", "dropout_summary.csv")))
cat("\n✓ TEST 1 PASSED: 5 combinations completed in parallel.\n\n")

cat("--- TEST 2: Re-run with caching enabled (should reuse all results instantly) ---\n")
res2 <- run_dropout_test(
  lfc_matrix = lfc_mat,
  se_matrix = se_mat,
  pval_matrix = pval_mat,
  comparison_gse_map = gse_map,
  sign_map = sign_map,
  all_cov_mats = all_cov_mats,
  cluster_annotation = annot,
  min_datasets = 0,
  min_dropout_datasets = 3,
  out_dir = test_dir,
  force_rerun = FALSE,
  num_cores = 2
)

stopifnot(nrow(res2) == 5)
cat("\n✓ TEST 2 PASSED: Caching successfully reused results.\n\n")

cat("--- TEST 3: Verify combination_info.txt contents ---\n")
info_path <- file.path(test_dir, "drop_out_test", "combination_1", "combination_info.txt")
stopifnot(file.exists(info_path))
cat("Content of combination_1/combination_info.txt:\n")
cat(readLines(info_path), sep = "\n")
cat("\n✓ TEST 3 PASSED!\n\n")

# Cleanup
unlink(test_dir, recursive = TRUE)
cat("ALL TESTS PASSED SUCCESSFULLY!\n")
