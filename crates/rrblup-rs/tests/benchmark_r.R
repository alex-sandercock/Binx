#!/usr/bin/env Rscript
# Benchmark R/rrBLUP mixed.solve and kin.blup

library(rrBLUP)

cat("R/rrBLUP Benchmark\n")
cat("==================\n\n")

# ============================================================
# mixed.solve with Z (marker matrix) - the common genomic case
# ============================================================
cat("mixed.solve with marker matrix Z (n genotypes x m markers):\n")
cat("------------------------------------------------------------\n")

benchmark_z <- function(name, expr, times = 3) {
  timings <- numeric(times)
  for (i in 1:times) {
    start <- Sys.time()
    result <- eval(expr)
    end <- Sys.time()
    timings[i] <- as.numeric(difftime(end, start, units = "secs")) * 1000
  }
  cat(sprintf("%s: %.2f ms (mean of %d runs)\n", name, mean(timings), times))
  invisible(mean(timings))
}

# Varying markers with fixed genotypes
n <- 500  # genotypes
for (m in c(1000, 5000, 10000, 20000)) {
  set.seed(42)
  # Marker matrix (centered, scaled)
  Z <- matrix(sample(c(-1, 0, 1), n * m, replace = TRUE), n, m)
  Z <- scale(Z, center = TRUE, scale = FALSE)

  # Phenotype
  beta <- rnorm(m) * 0.01  # marker effects
  y <- Z %*% beta + rnorm(n) * 0.5

  benchmark_z(sprintf("n=%d, m=%d markers", n, m), quote({
    mixed.solve(y, Z = Z)
  }))
}

cat("\n")

# Varying genotypes with fixed markers
m <- 10000  # markers
for (n in c(100, 200, 500, 1000)) {
  set.seed(42)
  Z <- matrix(sample(c(-1, 0, 1), n * m, replace = TRUE), n, m)
  Z <- scale(Z, center = TRUE, scale = FALSE)

  beta <- rnorm(m) * 0.01
  y <- Z %*% beta + rnorm(n) * 0.5

  benchmark_z(sprintf("n=%d genotypes, m=%d", n, m), quote({
    mixed.solve(y, Z = Z)
  }))
}

cat("\n")

# Benchmark function
benchmark <- function(name, expr, times = 5) {
  timings <- numeric(times)
  for (i in 1:times) {
    start <- Sys.time()
    result <- eval(expr)
    end <- Sys.time()
    timings[i] <- as.numeric(difftime(end, start, units = "secs")) * 1000
  }
  cat(sprintf("%s: %.2f ms (mean of %d runs, range: %.2f - %.2f ms)\n", 
              name, mean(timings), times, min(timings), max(timings)))
  invisible(mean(timings))
}

# (K-based benchmarks skipped for speed - see earlier results)

# ============================================================
# kin.blup benchmarks (reduced set)
# ============================================================
cat("kin.blup benchmarks:\n")
cat("--------------------\n")

for (n_geno in c(100, 200)) {
  set.seed(42)
  n_obs <- n_geno * 3  # 3 observations per genotype on average
  
  # Generate kinship matrix
  G <- matrix(rnorm(n_geno * 30), n_geno, 30)
  K <- tcrossprod(G) / 30
  K <- K + diag(n_geno) * 0.1
  geno_names <- paste0("G", 1:n_geno)
  rownames(K) <- colnames(K) <- geno_names
  
  # Generate phenotypes
  geno_effects <- rnorm(n_geno)
  obs_geno <- sample(geno_names, n_obs, replace = TRUE)
  y <- 5 + geno_effects[match(obs_geno, geno_names)] * 0.5 + rnorm(n_obs) * 0.3
  
  data <- data.frame(gid = obs_geno, y = y)
  
  benchmark(sprintf("n_geno=%d, n_obs=%d", n_geno, n_obs), quote({
    kin.blup(data, geno = "gid", pheno = "y", K = K)
  }), times = 3)
}

cat("\nDone.\n")
