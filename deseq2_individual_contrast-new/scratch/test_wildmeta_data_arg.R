# Test wildmeta with data argument approach
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

dat <- data.frame(
  yi = yi_vals,
  vi = vi_vals,
  gse = gse,
  comp = comp
)

# Test 1: with data= argument (which wildmeta internally expects)
cat("=== Test 1: rma.mv with data= argument ===\n")
m1 <- rma.mv(yi = yi, V = vi, random = ~ 1 | gse / comp, data = dat, method = "REML")
cat("Model fitted OK\n")
cat("call$data:", deparse(m1$call$data), "\n")

cwb1 <- tryCatch({
  Wald_test_cwb(m1, constrain_zero(1), R = 49, adjust = "CR2", seed = 1)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(cwb1)) {
  cat("SUCCESS!\n")
  print(cwb1)
}

# Test 2: with V as a matrix and data= argument
cat("\n=== Test 2: rma.mv with V as matrix + data= ===\n")
V_mat <- vcalc(vi = vi_vals, cluster = gse, obs = comp, rho = 0.5)
# V must be passed directly (not from data frame column)
m2 <- rma.mv(yi = yi, V = V_mat, random = ~ 1 | gse / comp, data = dat, method = "REML")
cat("Model fitted OK\n")

cwb2 <- tryCatch({
  Wald_test_cwb(m2, constrain_zero(1), R = 49, adjust = "CR2", seed = 1)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(cwb2)) {
  cat("SUCCESS!\n")
  print(cwb2)
}

# Test 3: RE Stouffer model with V=Cor_mat, W=neff, data=
cat("\n=== Test 3: RE Stouffer model (V=Cor, W=neff) + data= ===\n")
Cor_mat <- diag(1, 8)
for (g in unique(gse)) {
  idx <- which(gse == g)
  if (length(idx) > 1) {
    for (r in idx) {
      for (cc in idx) {
        if (r != cc) Cor_mat[r, cc] <- 0.5
      }
    }
  }
}

neff <- c(2.5, 3, 2, 4, 3.5, 2, 2.5, 3)
dat$neff <- neff

m3 <- rma.mv(yi = yi, V = Cor_mat, W = neff, random = ~ 1 | gse / comp,
             data = dat, method = "REML")
cat("Model fitted OK\n")

cwb3 <- tryCatch({
  Wald_test_cwb(m3, constrain_zero(1), R = 49, adjust = "CR2", seed = 1)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(cwb3)) {
  cat("SUCCESS!\n")
  print(cwb3)
}

# Test 4: Full RE Stouffer structure (5 GSE, 14 comparisons)
cat("\n=== Test 4: Full RE Stouffer (5 GSE, 14 comparisons) + data= ===\n")
gse_full <- c(rep("GSE001", 3), rep("GSE002", 2), rep("GSE003", 4),
              rep("GSE004", 3), rep("GSE005", 2))
comp_full <- paste0("comp_", seq_along(gse_full))
neff_full <- c(2.5, 3, 2, 4, 3.5, 2, 2.5, 3, 2, 3.5, 4, 3, 2.5, 3)

p_vals <- c(0.01, 0.03, 0.15, 0.005, 0.02, 0.08, 0.04, 0.12, 0.2,
            0.001, 0.01, 0.05, 0.03, 0.06)
z_vals <- qnorm(1 - p_vals)

Cor_full <- diag(1, 14)
for (g in unique(gse_full)) {
  idx <- which(gse_full == g)
  if (length(idx) > 1) {
    for (r in idx) {
      for (cc in idx) {
        if (r != cc) Cor_full[r, cc] <- 0.5
      }
    }
  }
}

dat_full <- data.frame(
  yi = z_vals,
  gse = gse_full,
  comp = comp_full,
  neff = neff_full
)

m4 <- rma.mv(yi = yi, V = Cor_full, W = neff, random = ~ 1 | gse / comp,
             data = dat_full, method = "REML")
cat("Model fitted OK\n")

cwb4 <- tryCatch({
  Wald_test_cwb(m4, constrain_zero(1), R = 99, adjust = "CR2", seed = 42)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(cwb4)) {
  cat("SUCCESS!\n")
  print(cwb4)
}

# Test 5: Edge case - only 2 GSE
cat("\n=== Test 5: Edge case - 2 GSE only ===\n")
idx2 <- gse_full %in% c("GSE001", "GSE002")
dat_2gse <- data.frame(
  yi = z_vals[idx2],
  gse = gse_full[idx2],
  comp = comp_full[idx2],
  neff = neff_full[idx2]
)
Cor_2gse <- Cor_full[idx2, idx2]

m5 <- tryCatch({
  rma.mv(yi = yi, V = Cor_2gse, W = neff, random = ~ 1 | gse / comp,
         data = dat_2gse, method = "REML")
}, error = function(e) {
  cat("  rma.mv ERROR:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(m5)) {
  cat("  Model fitted OK\n")
  cwb5 <- tryCatch({
    Wald_test_cwb(m5, constrain_zero(1), R = 49, adjust = "CR2", seed = 42)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(cwb5)) {
    cat("  SUCCESS!\n")
    print(cwb5)
  }
}

cat("\n=== DONE ===\n")
