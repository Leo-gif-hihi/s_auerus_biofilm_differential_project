# ---------------------------------------------------------
# 1. SIMULATION UNDER THE NULL HYPOTHESIS (No True DEGs)
# ---------------------------------------------------------
n_genes <- 100000

# Under the Null Hypothesis, p-values are completely random.
# Mathematically, this means they follow a Uniform(0,1) distribution.
p_null <- runif(n_genes, min = 0, max = 1)

# Apply the two-sided Z-score transformation
z_null <- qnorm(1 - p_null)

# Plot the Null Hypothesis
png("Simulated_Null_Hypothesis.png", width = 800, height = 400, res = 120)
hist(z_null,
    breaks = 100, freq = FALSE,
    col = "lightgreen", border = "white",
    main = "Null Hypothesis: Matches N(0,1) Perfectly",
    xlab = "Z-score"
)

# Overlay theoretical N(0,1) curve
x_range <- seq(-5, 5, length.out = 500)
lines(x_range, dnorm(x_range), col = "red", lwd = 2)
dev.off()


# ---------------------------------------------------------
# 2. SIMULATION WITH BIOLOGICAL SIGNAL (Your Data)
# ---------------------------------------------------------
# Let's say 20% of genes are True DEGs (either up or down regulated)
n_true_degs <- 20000
n_null_genes <- 80000

# True DEGs have very small p-values (skewed towards 0)
p_degs <- rbeta(n_true_degs, shape1 = 0.5, shape2 = 10)

# Combine them with the random null genes
p_real <- c(runif(n_null_genes), p_degs)

# In your real data, you also cap the lowest p-values at the machine limit
p_real[p_real < .Machine$double.eps] <- .Machine$double.eps

# Apply the two-sided Z-score transformation
z_real <- qnorm(1 - p_real)

# Plot the Biological Signal
png("Simulated_Biological_Signal.png", width = 800, height = 400, res = 120)
hist(z_real,
    breaks = 100, freq = FALSE,
    col = "lightblue", border = "white",
    main = "Biological Signal: Matches Your Plot",
    xlab = "Z-score", xlim = c(-4, 9)
)

# Overlay theoretical N(0,1) curve
lines(x_range, dnorm(x_range), col = "red", lwd = 2)
dev.off()
